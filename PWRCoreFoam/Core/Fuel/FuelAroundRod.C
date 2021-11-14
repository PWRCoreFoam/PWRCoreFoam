/*----------------------------------------------------------------------------*\
| Class:                                                                       |
|     Foam::Fuel                                                               |
\*----------------------------------------------------------------------------*/


#include "Fuel.H"

#include "fvCFD.H"
#include "rhoThermo.H"
#include "rodFieldsFwd.H"
#include "subChannelMesh.H"


//****************************************************************************//

/*----------------------------------------------------------------------------*\
|                       Member functions                                       |
\*----------------------------------------------------------------------------*/

// Private member functions

scalar
Foam::Fuel::get_avg_T_around_rod(label edge_i)
{
  scalar retval = 0.0;

  forAll(mesh_.edgeCells()[edge_i], j)
  {
    label cell_j = mesh_.edgeCells()[edge_i][j];
    retval += thermo_.T()[cell_j];
  }

  retval = retval / mesh_.edgeCells()[edge_i].size();

  return retval;
}


scalar
Foam::Fuel::get_avg_U_z_around_rod(label edge_i)
{
  scalar retval = 0.0;

  forAll(mesh_.edgeCells()[edge_i], j)
  {
    label cell_j = mesh_.edgeCells()[edge_i][j];
    retval += U_[cell_j].z();
  }

  retval = retval / mesh_.edgeCells()[edge_i].size();

  return retval;
}


scalar
Foam::Fuel::get_avg_rho_around_rod(label edge_i)
{
  scalar retval = 0.0;

  forAll(mesh_.edgeCells()[edge_i], j)
  {
    label cell_j = mesh_.edgeCells()[edge_i][j];
    retval += thermo_.rho()()[cell_j];
  }

  retval = retval / mesh_.edgeCells()[edge_i].size();

  return retval;
}


scalar
Foam::Fuel::get_avg_Cp_around_rod(label edge_i)
{
  scalar retval = 0.0;

  forAll(mesh_.edgeCells()[edge_i], j)
  {
    label cell_j = mesh_.edgeCells()[edge_i][j];
    retval += Cp_[cell_j];
  }

  retval = retval / mesh_.edgeCells()[edge_i].size();

  return retval;
}


scalar
Foam::Fuel::get_avg_mu_around_rod(label edge_i)
{
  scalar retval = 0.0;

  forAll(mesh_.edgeCells()[edge_i], j)
  {
    label cell_j = mesh_.edgeCells()[edge_i][j];
    retval += thermo_.mu()[cell_j];
  }

  retval = retval / mesh_.edgeCells()[edge_i].size();

  return retval;
}


scalar
Foam::Fuel::get_avg_kappa_around_rod(label edge_i)
{
  scalar retval = 0.0;

  forAll(mesh_.edgeCells()[edge_i], j)
  {
    label cell_j = mesh_.edgeCells()[edge_i][j];
    retval += kappa_[cell_j];
  }

  retval = retval / mesh_.edgeCells()[edge_i].size();

  return retval;
}


scalar
Foam::Fuel::get_avg_alpha_around_rod(label edge_i)
{
  scalar retval = 0.0;

  forAll(mesh_.edgeCells()[edge_i], j)
  {
    label cell_j = mesh_.edgeCells()[edge_i][j];
    retval += thermo_.alpha()[cell_j];
  }

  retval = retval / mesh_.edgeCells()[edge_i].size();

  return retval;
}


