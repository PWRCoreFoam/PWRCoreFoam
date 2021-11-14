/*----------------------------------------------------------------------------*\
| Class:                                                                       |
|     Foam::GeometricRodField                                                  |
\*----------------------------------------------------------------------------*/

#include "data.H"
#include "demandDrivenData.H"
#include "dictionary.H"
#include "edgeList.H"
#include "GeometricRodField.H"
#include "OFstream.H"
#include "pointField.H"
#include "Time.H"


//****************************************************************************//

/*----------------------------------------------------------------------------*\
|                       Constructors                                           |
\*----------------------------------------------------------------------------*/

template<class Type, class GeoMesh>
Foam::GeometricRodField<Type, GeoMesh>::GeometricRodField
(
    const IOobject& io,
    const Mesh& mesh,
    const dimensioned<Type>& dt
)
:
    DimensionedField<Type, GeoMesh>(io, mesh, dt, false),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(NULL),
    fieldPrevIterPtr_(NULL),
    subChannelMesh_(mesh),
    rodEdgesList_(mesh.get_rodEdgesList()),
    rodEdgesListOnAssemblyEdge_(mesh.get_rodEdgesListOnAssemblyEdge()),
    rodEdgesListOnCoreEdge_(mesh.get_rodEdgesListOnCoreEdge()),
    rodEdgesLengthList_(mesh.get_rodEdgesLengthList()),
    rodEdgesCenterList_(mesh.get_rodEdgesCenterList()),
    z_min(GREAT),
    z_max(-1.0*GREAT)
{

    if (debug)
    {
        Info<< "GeometricField<Type, PatchField, GeoMesh>::GeometricField : "
               "creating temporary"
            << endl << this->info() << endl;
    }

    calc_z_min_max();
}


template<class Type, class GeoMesh>
Foam::GeometricRodField<Type, GeoMesh>::GeometricRodField
(
    const IOobject& io,
    const Mesh& mesh
)
:
    DimensionedField<Type, GeoMesh>(io, mesh, dimless, false),
    timeIndex_(this->time().timeIndex()),
    field0Ptr_(NULL),
    fieldPrevIterPtr_(NULL),
    subChannelMesh_(mesh),
    rodEdgesList_(mesh.get_rodEdgesList()),
    rodEdgesListOnAssemblyEdge_(mesh.get_rodEdgesListOnAssemblyEdge()),
    rodEdgesListOnCoreEdge_(mesh.get_rodEdgesListOnCoreEdge()),
    rodEdgesLengthList_(mesh.get_rodEdgesLengthList()),
    rodEdgesCenterList_(mesh.get_rodEdgesCenterList()),
    z_min(GREAT),
    z_max(-1.0*GREAT)
{

    if (debug)
    {
        Info<< "Finishing read-construct of "
               "GeometricField<Type, PatchField, GeoMesh>"
            << endl << this->info() << endl;
    }

    calc_z_min_max();
}


/*----------------------------------------------------------------------------*\
|                       Destructor                                             |
\*----------------------------------------------------------------------------*/

template<class Type, class GeoMesh>
Foam::GeometricRodField<Type, GeoMesh>::~GeometricRodField()
{
    deleteDemandDrivenData(field0Ptr_);
    deleteDemandDrivenData(fieldPrevIterPtr_);
}


/*----------------------------------------------------------------------------*\
|                       Member functions                                       |
\*----------------------------------------------------------------------------*/

template<class Type, class GeoMesh>
void Foam::GeometricRodField<Type, GeoMesh>::calc_z_min_max()
{
  forAll(rodEdgesList_, i)
  {
    label edge_i = rodEdgesList_[i];

    label startLabel = subChannelMesh_.edges()[edge_i].start();
    vector startPoint = subChannelMesh_.points()[startLabel];

    z_min = (startPoint.z() < z_min) ? startPoint.z() : z_min;
    z_max = (startPoint.z() > z_max) ? startPoint.z() : z_max;

    label endLabel = subChannelMesh_.edges()[edge_i].end();
    vector endPoint = subChannelMesh_.points()[endLabel];

    z_min = (endPoint.z() < z_min) ? endPoint.z() : z_min;
    z_max = (endPoint.z() > z_max) ? endPoint.z() : z_max;
  }
}


template<class Type, class GeoMesh>
void Foam::GeometricRodField<Type, GeoMesh>::writeRodField() const
{

    OFstream os((this->time().timePath())/(this->name())+".vtk");

    // Write title
    os << "# vtk DataFile Version 2.0" << endl
       << "Rod Data"                   << endl
       << "ASCII"                      << endl
       << "DATASET UNSTRUCTURED_GRID"  << endl;


    label rod_cells_number = rodEdgesList_.size();

    scalar pitch_ = rodMesh().get_pitch();
    scalar rodDelta(pitch_ / 8.0);

    // Write points
    os << "POINTS        " << rod_cells_number*8 << " double" << endl;

    forAll(rodEdgesList_, i)
    {
      label edge_i = rodEdgesList_[i];

      label p1_label = subChannelMesh_.edges()[edge_i].start();
      vector p1 = subChannelMesh_.points()[p1_label];

      label p2_label = subChannelMesh_.edges()[edge_i].end();
      vector p2 = subChannelMesh_.points()[p2_label];

      if (p1.z() > p2.z())
      {
        Foam::Swap(p1, p2);
      }

      os << p1.x() - rodDelta << "    "
         << p1.y() - rodDelta << "    "
         << p1.z() << endl;

      os << p1.x() + rodDelta << "    "
         << p1.y() - rodDelta << "    "
         << p1.z() << endl;

      os << p1.x() - rodDelta << "    "
         << p1.y() + rodDelta << "    "
         << p1.z() << endl;

      os << p1.x() + rodDelta << "    "
         << p1.y() + rodDelta << "    "
         << p1.z() << endl;

      os << p2.x() - rodDelta << "    "
         << p2.y() - rodDelta << "    "
         << p2.z() << endl;

      os << p2.x() + rodDelta << "    "
         << p2.y() - rodDelta << "    "
         << p2.z() << endl;

      os << p2.x() - rodDelta << "    "
         << p2.y() + rodDelta << "    "
         << p2.z() << endl;

      os << p2.x() + rodDelta << "    "
         << p2.y() + rodDelta << "    "
         << p2.z() << endl;
    }

    // Write cells
    os << "CELLS    " << rod_cells_number << " " << 9*rod_cells_number << endl;

    forAll(rodEdgesList_, i)
    {
      os << "   8   " << i*8+0 << "      " << i*8+1 << "      "
                      << i*8+2 << "      " << i*8+3 << "      "
                      << i*8+4 << "      " << i*8+5 << "      "
                      << i*8+6 << "      " << i*8+7 << "      " << endl;
    }

    os << "CELL_TYPES    " << rod_cells_number << endl;

    for (label i=0; i<rod_cells_number; i++)
    {
      os << 11 << endl;
    }

    // Write cells data
    os << "CELL_DATA    " << rod_cells_number << endl;
    os << "SCALARS    " << this->name() << "  float" << endl;
    os << "LOOKUP_TABLE default" << endl;

    forAll(rodEdgesList_, i)
    {
      os << this->operator[](i) << endl;
    }
}


template<class Type, class GeoMesh>
void Foam::GeometricRodField<Type, GeoMesh>::writeRodOnAssemblyEdge() const
{

    OFstream os((this->time().timePath())/(this->name())+"_Rods_on_Assembly.vtk");

    // Write title
    os << "# vtk DataFile Version 2.0" << endl
       << "Rod Data"                   << endl
       << "ASCII"                      << endl
       << "DATASET UNSTRUCTURED_GRID"  << endl;


    label rod_cells_number = rodEdgesListOnAssemblyEdge_.size();

    scalar pitch_ = rodMesh().get_pitch();
    scalar rodDelta(pitch_ / 12.0);

    // Write points
    os << "POINTS        " << rod_cells_number*8 << " double" << endl;

    forAll(rodEdgesListOnAssemblyEdge_, i)
    {
      label edge_i = rodEdgesListOnAssemblyEdge_[i];

      label p1_label = subChannelMesh_.edges()[edge_i].start();
      vector p1 = subChannelMesh_.points()[p1_label];

      label p2_label = subChannelMesh_.edges()[edge_i].end();
      vector p2 = subChannelMesh_.points()[p2_label];

      if (p1.z() > p2.z())
      {
        Foam::Swap(p1, p2);
      }

      os << p1.x() - rodDelta << "    "
         << p1.y() - rodDelta << "    "
         << p1.z() << endl;

      os << p1.x() + rodDelta << "    "
         << p1.y() - rodDelta << "    "
         << p1.z() << endl;

      os << p1.x() - rodDelta << "    "
         << p1.y() + rodDelta << "    "
         << p1.z() << endl;

      os << p1.x() + rodDelta << "    "
         << p1.y() + rodDelta << "    "
         << p1.z() << endl;

      os << p2.x() - rodDelta << "    "
         << p2.y() - rodDelta << "    "
         << p2.z() << endl;

      os << p2.x() + rodDelta << "    "
         << p2.y() - rodDelta << "    "
         << p2.z() << endl;

      os << p2.x() - rodDelta << "    "
         << p2.y() + rodDelta << "    "
         << p2.z() << endl;

      os << p2.x() + rodDelta << "    "
         << p2.y() + rodDelta << "    "
         << p2.z() << endl;
    }

    // Write cells
    os << "CELLS    " << rod_cells_number << " " << 9*rod_cells_number << endl;

    forAll(rodEdgesListOnAssemblyEdge_, i)
    {
      os << "   8   " << i*8+0 << "      " << i*8+1 << "      "
                      << i*8+2 << "      " << i*8+3 << "      "
                      << i*8+4 << "      " << i*8+5 << "      "
                      << i*8+6 << "      " << i*8+7 << "      " << endl;
    }

    os << "CELL_TYPES    " << rod_cells_number << endl;

    for (label i=0; i<rod_cells_number; i++)
    {
      os << 11 << endl;
    }

    // Write cells data
    os << "CELL_DATA    " << rod_cells_number << endl;
    os << "SCALARS    " << this->name() << "  float" << endl;
    os << "LOOKUP_TABLE default" << endl;

    forAll(rodEdgesListOnAssemblyEdge_, i)
    {
      os << scalar(1.0) << endl;
    }
}


template<class Type, class GeoMesh>
void Foam::GeometricRodField<Type, GeoMesh>::writeRodOnCoreEdge() const
{

    OFstream os((this->time().timePath())/(this->name())+"_Rods_on_Core.vtk");

    // Write title
    os << "# vtk DataFile Version 2.0" << endl
       << "Rod Data"                   << endl
       << "ASCII"                      << endl
       << "DATASET UNSTRUCTURED_GRID"  << endl;


    label rod_cells_number = rodEdgesListOnCoreEdge_.size();

    scalar pitch_ = rodMesh().get_pitch();
    scalar rodDelta(pitch_ / 8.0);

    // Write points
    os << "POINTS        " << rod_cells_number*8 << " double" << endl;

    forAll(rodEdgesListOnCoreEdge_, i)
    {
      label edge_i = rodEdgesListOnCoreEdge_[i];

      label p1_label = subChannelMesh_.edges()[edge_i].start();
      vector p1 = subChannelMesh_.points()[p1_label];

      label p2_label = subChannelMesh_.edges()[edge_i].end();
      vector p2 = subChannelMesh_.points()[p2_label];

      if (p1.z() > p2.z())
      {
        Foam::Swap(p1, p2);
      }

      os << p1.x() - rodDelta << "    "
         << p1.y() - rodDelta << "    "
         << p1.z() << endl;

      os << p1.x() + rodDelta << "    "
         << p1.y() - rodDelta << "    "
         << p1.z() << endl;

      os << p1.x() - rodDelta << "    "
         << p1.y() + rodDelta << "    "
         << p1.z() << endl;

      os << p1.x() + rodDelta << "    "
         << p1.y() + rodDelta << "    "
         << p1.z() << endl;

      os << p2.x() - rodDelta << "    "
         << p2.y() - rodDelta << "    "
         << p2.z() << endl;

      os << p2.x() + rodDelta << "    "
         << p2.y() - rodDelta << "    "
         << p2.z() << endl;

      os << p2.x() - rodDelta << "    "
         << p2.y() + rodDelta << "    "
         << p2.z() << endl;

      os << p2.x() + rodDelta << "    "
         << p2.y() + rodDelta << "    "
         << p2.z() << endl;
    }

    // Write cells
    os << "CELLS    " << rod_cells_number << " " << 9*rod_cells_number << endl;

    forAll(rodEdgesListOnCoreEdge_, i)
    {
      os << "   8   " << i*8+0 << "      " << i*8+1 << "      "
                      << i*8+2 << "      " << i*8+3 << "      "
                      << i*8+4 << "      " << i*8+5 << "      "
                      << i*8+6 << "      " << i*8+7 << "      " << endl;
    }

    os << "CELL_TYPES    " << rod_cells_number << endl;

    for (label i=0; i<rod_cells_number; i++)
    {
      os << 11 << endl;
    }

    // Write cells data
    os << "CELL_DATA    " << rod_cells_number << endl;
    os << "SCALARS    " << this->name() << "  float" << endl;
    os << "LOOKUP_TABLE default" << endl;

    forAll(rodEdgesListOnCoreEdge_, i)
    {
      os << scalar(2.0) << endl;
    }
}


template<class Type, class GeoMesh>
void Foam::GeometricRodField<Type, GeoMesh>::writeRodField_Test() const
{

    OFstream os((this->time().timePath())/(this->name())+"_Test.vtk");

    os << "# vtk DataFile Version 2.0" << endl
       << "Rod Data"                   << endl
       << "ASCII"                      << endl
       << "DATASET UNSTRUCTURED_GRID"  << endl;

    labelList rods = rodMesh().get_rodEdgesList();

    label nRodCells = rods.size();

    scalar pitch_ = rodMesh().get_pitch();
    scalar rodDelta(pitch_ / 8.0);

    const pointField & all_points = rodMesh().points();
    const edgeList & all_edges = rodMesh().edges();

    vector p1;
    vector p2;

    scalar z_max(-1.0*GREAT);
    scalar z_min(GREAT);


    os << "POINTS        " << nRodCells*8 << " double" << endl;

    for (label i=0; i<nRodCells; i++)
    {
      label edgeI = rods[i];

      p1 = all_points[all_edges[edgeI][0]];
      p2 = all_points[all_edges[edgeI][1]];

      if (p1.z() > p2.z())
      {
        Foam::Swap(p1, p2);
      }

      os << p1.x() - rodDelta << "    "
         << p1.y() - rodDelta << "    "
         << p1.z() << endl;

      os << p1.x() + rodDelta << "    "
         << p1.y() - rodDelta << "    "
         << p1.z() << endl;

      os << p1.x() - rodDelta << "    "
         << p1.y() + rodDelta << "    "
         << p1.z() << endl;

      os << p1.x() + rodDelta << "    "
         << p1.y() + rodDelta << "    "
         << p1.z() << endl;

      if (p1.z() > z_max)
      {
        z_max = p1.z();
      }
      if (p1.z() < z_min)
      {
        z_min = p1.z();
      }

      os << p2.x() - rodDelta << "    "
         << p2.y() - rodDelta << "    "
         << p2.z() << endl;

      os << p2.x() + rodDelta << "    "
         << p2.y() - rodDelta << "    "
         << p2.z() << endl;

      os << p2.x() - rodDelta << "    "
         << p2.y() + rodDelta << "    "
         << p2.z() << endl;

      os << p2.x() + rodDelta << "    "
         << p2.y() + rodDelta << "    "
         << p2.z() << endl;

      if (p2.z() > z_max)
      {
        z_max = p2.z();
      }
      if (p2.z() < z_min)
      {
        z_min = p2.z();
      }


    }

    os << "CELLS    " << nRodCells << " " << 9*nRodCells << endl;

    for (label i=0; i<nRodCells; i++)
    {
      os << "   8   " << i*8+0 << "      "
                      << i*8+1 << "      "
                      << i*8+2 << "      "
                      << i*8+3 << "      "
                      << i*8+4 << "      "
                      << i*8+5 << "      "
                      << i*8+6 << "      "
                      << i*8+7 << "      " << endl;
    }

    os << "CELL_TYPES    " << nRodCells << endl;

    for (label i=0; i<nRodCells; i++)
    {
      os << 11 << endl;
    }

    os << "CELL_DATA    " << nRodCells << endl;

    os << "SCALARS    " << this->name() << "  float" << endl;
    os << "LOOKUP_TABLE default" << endl;


    for (label i=0; i<nRodCells; i++)
    {
      label edgeI = rods[i];

      p1 = all_points[all_edges[edgeI][0]];
      p2 = all_points[all_edges[edgeI][1]];

      scalar z_((p1.z()+p2.z())/2.0);

      os << sin((z_max - z_)/(z_max - z_min)*3.1415926) << endl;

    }

}


//****************************************************************************//


