#include "distance_marcher.h"
#include "math.h"
#include <stdexcept>
#include <map>
#include <set>
//TODO: Remove iostream
#include <iostream>

double distanceMarcher::updatePointOrderOne(int i)
{
  double a,b,c;
  a=b=c=0;
  int naddr=0;
  for (int dim=0; dim<dim_; dim++)
  {
    double value = maxDouble;
    for (int j=-1; j<2; j+=2) // each direction
    {
      naddr = _getN(i,dim,j,Mask);
      if (naddr!=-1 && flag_[naddr]==Frozen)
      {
        if (fabs(distance_[naddr])<fabs(value))
        {
          value = distance_[naddr];
        }
      }
    }
    if (value<maxDouble)
    {
      a+=idx2_[dim];
      b-=idx2_[dim]*2*value;
      c+=idx2_[dim]*pow(value,2);
    }
  }
  try {
    return solveQuadratic(i,a,b,c);
  } catch(std::runtime_error& err) {
    // if the quadratic equation can't be solved, use the
    // position of the minimum as a more reasonable approximation
    return -b / (2.0 * a);
  }
}


// second order point update
// update the distance from the frozen points
const double aa         =  9.0/4.0;
const double oneThird   =  1.0/3.0;
double distanceMarcher::updatePointOrderTwo(int i)
{
  double a,b,c;
  a=b=c=0;
  int naddr=0;
  for (int dim=0; dim<dim_; dim++)
  {
    double value1 = maxDouble;
    double value2 = maxDouble;
    for (int j=-1; j<2; j+=2) // each direction
    {
      naddr = _getN(i,dim,j,Mask);
      if (naddr!=-1 && flag_[naddr]==Frozen)
      {
        if (fabs(distance_[naddr])<fabs(value1))
        {
          value1 = distance_[naddr];
          int naddr2 = _getN(i,dim,j*2,Mask);
          if (naddr2!=-1 &&
              flag_[naddr2]==Frozen &&
              ((distance_[naddr2]<=value1 && value1 >=0) ||
               (distance_[naddr2]>=value1 && value1 <=0)))
          {
            value2=distance_[naddr2];
          }
        }
      }
    }
    if (value2<maxDouble)
    {
      double tp = oneThird*(4*value1-value2);
      a+=idx2_[dim]*aa;
      b-=idx2_[dim]*2*aa*tp;
      c+=idx2_[dim]*aa*pow(tp,2);
    }
    else if (value1<maxDouble)
    {
      a+=idx2_[dim];
      b-=idx2_[dim]*2*value1;
      c+=idx2_[dim]*pow(value1,2);
    }
  }
  try {
    return solveQuadratic(i,a,b,c);
  } catch(std::runtime_error& err) {
    // if the second order method fails, try the first order method instead
    return updatePointOrderOne(i);
  }
}


// find and return the correct root
double distanceMarcher::solveQuadratic(int i, const double &a,
                                       const double &b,
                                       double &c)
{
  c-=1;
  double det = pow(b,2)-4*a*c;
  if (det >= 0)
  {
    if (phi_[i] > doubleEpsilon) { return (-b + sqrt(det)) / 2.0 / a; }
    else                         { return (-b - sqrt(det)) / 2.0 / a; }
  }
  else
  {
    throw std::runtime_error("Negative discriminant in distance marcher quadratic.");
  }
}


void distanceMarcher::initalizeFrozen2()
{
  //loop over phi to find zero values
  //  and mark them as frozen
  for (int i=0; i<size_; i++)
  {
    if (flag_[i] != Mask && phi_[i]==0.0)
    {
      flag_[i]=Frozen;
      distance_[i]=0.0;
    }
  }
  //loop over all of phi and for each point check each direction
  //  to see if we cross the zero level set
  //     if so calculate the minimum distance to the zero level set
  //     mark as frozen.
  for (int i=0; i<size_; i++)
  if (flag_[i] == Far)
  {
    double ldistance[MaximumDimension];
    bool borders=false;
    for (int dim=0; dim<dim_; dim++)
    {
      ldistance[dim]=0;
      for (int j=-1; j<2; j+=2) // each direction
      {
        int naddr = _getN(i,dim,j,Mask);
        if (naddr!=-1 && phi_[i] * phi_[naddr]<0)
        {
          // this cell and neighbor span the zero level set.
          borders=true;
          //calculate the distance to the zero level set.
          double d = dx_[dim]*phi_[i]/(phi_[i]-phi_[naddr]);
          if (ldistance[dim]==0 || ldistance[dim]>d)
          {
            ldistance[dim] = d;
          }
        }
      } // for each direction
    } // for each dimension
    if (borders)
    {
      double dsum = 0;
      for (int dim=0; dim<dim_; dim++)
        if (ldistance[dim]>0) dsum += 1/ldistance[dim]/ldistance[dim];
      if (phi_[i]<0)
        distance_[i] = -sqrt(1/dsum);
      else distance_[i] = sqrt(1/dsum);
      flag_[i]=Frozen;
    }
  }// for each point in the far field
}

void distanceMarcher::initalizeFrozen()
{
  //loop over phi to find zero values
  //  and mark them as frozen
  for (int i=0; i<size_; i++)
  {
    if (flag_[i] != Mask && phi_[i]==0.0)
    {
      flag_[i]=Frozen;
      distance_[i]=0.0;
    }
  }
  // loop over all of phi and for each point check each direction to see if we
  // cross the zero level set. If we do, then the ratio between phi at these
  // points must be fixed so that the zero level-set is anchored.
  //
  // Looping over each edge between nodes i and j (with i < j) which crosses the zero-level-set,
  // add an edge in a multimap with key i, value (j, r) where r is the ratio between the value at j and at i.
  // If there is an edge between i and j and one between j and k, then the ratio between i and k is r_ij * r_jk.
  //
  // When all edges are added, you can pick any node and compute what the
  // relative sizes of all nodes that directly or indirectly depend on this
  // node by traversing the edges. Now we know what the relative sizes of some
  // chunk of points should be. How do we select with which factor we should
  // multiply this? Let's pick the factor that minimizes the error. If the
  // relative sizes is stored in a vector r and the actual distances in a
  // vector d, then the closest approximation of d along r is (d . r / |r|^2) * r.

  // DestinationRatio
  typedef std::pair<int, double> DestRat;

  // Forward / backward maps
  std::multimap<int, DestRat> fw_map;
  std::multimap<int, DestRat> bw_map;

  // Distance map
  std::map<int, double> d_map;

  // Loop over non-frozen values
  for (int i=0; i<size_; i++)
    if (flag_[i] == Far)
    {
      for (int dim=0; dim<dim_; dim++)
      {
        // Since we only need to check each edge once (i.e. not i-j and j-i)
        // we can just check the positive direction
        int j = _getN(i,dim,+1,Mask);

        if (j !=-1 && phi_[i] * phi_[j]<0)
        {
          // this cell and neighbor span the zero level set.

          //calculate the distance to the zero level set.
          double d = dx_[dim]*phi_[i]/(phi_[i]-phi_[j]);

          // Calculate the ratio between the distance values of node i and j
          double r = phi_[j]/phi_[i];
          double r_inv = phi_[i]/phi_[j];

          // Insert edge from both i->j and j->i in both the backward map and forward map.
          fw_map.insert({i, DestRat(j, r)});
          bw_map.insert({j, DestRat(i, r_inv)});

          // Update the distance map if needed.
          auto d_it = d_map.find(i);
          if (d_it == d_map.end() || abs(d_it->second) > d) {
            if (phi_[i] < 0)
              d_map[i]= -d;
            else
              d_map[i] = d;
          }
          d_it = d_map.find(j);
          if (d_it == d_map.end() || abs(d_it->second) > dx_[dim]-d) {
            if (phi_[j] < 0)
              d_map[j] = -dx_[dim]+d;
            else
              d_map[j] = dx_[dim]-d;
          }

          // Mark both i and j as "to be frozen"
          flag_[i]=Frozen;
          flag_[j]=Frozen;
        }
      } // for each dimension
    }// for each point in the far field

  // Now we go through the nodes and set the distances
  // For that we need to keep track of which nodes we've already updated.
  std::set<int> visited;

  // Loop over all where there are distances noted
  // (i.e. all the ones that were frozen in the previous loop)
  for (auto it = d_map.begin(); it != d_map.end(); ++it) {

    // Pick i as the "root" node
    // The relations to all nodes dependant on this one will be computed
    // relative to this one. If this node was already visited, it is
    // because it was dependent on some other node, thus we skip it
    int i = it->first;
    if (visited.find(i) != visited.end())
      continue;

    // Start with a relative factor of 1
    double r = 1;

    // Create a map that stores the nodes that are dependent on i but we have not visited,
    // as well as their relative factors.
    std::map<int, double> to_visit;
    // Create a map that stores the nodes that are dependent on i and their relative factors,
    // for all the visited nodes.
    std::map<int, double> null_space_basis_vector;

    // First entry, i and 1.
    visited.insert(i);
    null_space_basis_vector.insert({i, r});

    // Add all of the nodes directly dependent on i as to-be-visited, with their relative factors
    auto fw_range = fw_map.equal_range(i);
    for (auto jt = fw_range.first; jt != fw_range.second; ++jt){
      DestRat dest_rat = jt->second;

      int j = dest_rat.first;

      if (visited.find(j) == visited.end())
        to_visit.insert({j, r * dest_rat.second});
    }
    auto bw_range = bw_map.equal_range(i);
    for (auto jt = bw_range.first; jt != bw_range.second; ++jt){
      DestRat dest_rat = jt->second;

      int j = dest_rat.first;

      if (visited.find(j) == visited.end())
        to_visit.insert({j, r * dest_rat.second});
    }

    // We continue like this until there are no more nodes to be visited
    while (!to_visit.empty()) {
      // Get any node from the to_visit store, extract index and relative factor
      auto next_it = to_visit.begin();
      int i = next_it->first;
      double r = next_it->second;

      // Mark as visited
      visited.insert(i);

      // Insert into basis-vector
      null_space_basis_vector.insert({i, r});

      // Add all nodes directly dependent on *this* node to the to_visit map
      auto fw_range = fw_map.equal_range(i);
      for (auto jt = fw_range.first; jt != fw_range.second; ++jt){
        DestRat dest_rat = jt->second;

        int j = dest_rat.first;

        if (visited.find(j) == visited.end())
          to_visit.insert({j, r * dest_rat.second});
      }
      auto bw_range = bw_map.equal_range(i);
      for (auto jt = bw_range.first; jt != bw_range.second; ++jt){
        DestRat dest_rat = jt->second;

        int j = dest_rat.first;

        if (visited.find(j) == visited.end())
          to_visit.insert({j, r * dest_rat.second});
      }

      // Erase the current entry and continue with the next
      to_visit.erase(next_it);

    } // while to visit isn't empty


    // Compute scalar product with d_map, and norm of null space basis vector
    double norm_sq = 0.0;
    double scal_prod = 0.0;
    for (auto jt = null_space_basis_vector.begin(); jt != null_space_basis_vector.end(); ++jt) {
      norm_sq += jt->second * jt->second;
      scal_prod += jt->second * d_map[jt->first];
    }
    // Set distances accordingly: d[i] = (nsbv . d_map) / |nsbv|^2 * nsbv
    for (auto jt = null_space_basis_vector.begin(); jt != null_space_basis_vector.end(); ++jt) {
      distance_[jt->first] = scal_prod / norm_sq * jt->second;

    }
  } // For all points in d_map
}
