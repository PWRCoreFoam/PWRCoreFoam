/*----------------------------------------------------------------------------*\
| Class:                                                                       |
|     Foam::subChannelMesh                                                     |
\*----------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "slicedSurfaceFields.H"
#include "slicedVolFields.H"
#include "subChannelMesh.H"


//****************************************************************************//

namespace Foam
{

defineTypeNameAndDebug(subChannelMesh, 0);

}


/*----------------------------------------------------------------------------*\
|                       Constructors                                           |
\*----------------------------------------------------------------------------*/

Foam::subChannelMesh::subChannelMesh(const IOobject & io)
:
  fvMesh(io),
  rodEdgesListPtr_(NULL),
  rodEdgesLengthListPtr_(NULL),
  rodEdgesCenterListPtr_(NULL),
  rodEdgesListOnAssemblyEdgePtr_(NULL),
  rodEdgesListOnCoreEdgePtr_(NULL),
  gapListPtr_(NULL),
  CoreMapDimension_(0),
  AssemblyRodsDimension_(0),
  AssemblyWidth_(0.0),
  CenterAssemblyXIndex_(-1),
  CenterAssemblyYIndex_(-1),
  center_rod_X_index_(-1.0),
  center_rod_Y_index_(-1.0),
  rod_assembly_index_X_Y_(NULL),
  rod_rod_index_X_Y_(NULL),
  core_type_(NULL),
  cellRodsNumberList_(NULL),
  cellRodsList_(NULL),
  RodOuterDiameter_(0.0),
  pitch_(0.0),
  correct_rod_arrangement_(Switch("no")),
  test_case_with_assembly_box_(Switch("no")),
  assembly_box_corner_radius_(0.0),
  mixing_coefficient_(0.0),
  subChannelMeshBuildOrNot(false),
  rodsBuildOrNot(false),
  mesh_modified_or_not(false),
  flow_area_(this->nCells(), 0.0),
  hydraulic_diameter_(this->nCells(), 0.0)
{
  readReactorCore();
}


/*----------------------------------------------------------------------------*\
|                       Destructor                                             |
\*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*\
|                       Member functions                                       |
\*----------------------------------------------------------------------------*/

void
Foam::subChannelMesh::modify_mesh()
{

  if (!rodsBuildOrNot)
  {
    FatalErrorIn("void Foam::subChannelMesh::modify_mesh()")
              << "rods should be built first !"
              << exit(FatalError);
  }

  labelList rodsSectionList_ = *rodEdgesListPtr_;


  /*--------------------------------------------------------------------------*\
  | cells                                                                      |
  \*--------------------------------------------------------------------------*/

  volScalarField::DimensionedInternalField * vptr = get_VPtr();
  if (!vptr)
  {
    FatalErrorIn("void Foam::subChannelMesh::modify_mesh()")
              << "Can't get cell volume field !"
              << exit(FatalError);
  }

  volScalarField::DimensionedInternalField & cell_volumes = *vptr;

  forAll(rodsSectionList_, i)
  {
    label edge_i = (*rodEdgesListPtr_)[i];
    labelList rodsCellList_ = edgeCells()[edge_i];

    scalar pi_ = 3.14159;
    scalar rod_outer_diameter_ = RodOuterDiameter_;
    scalar rod_section_length_ = get_rod_section_length(edge_i);
    scalar rod_section_volume_ = 1.0 / 4.0 * pi_ * ::sqr(rod_outer_diameter_) * rod_section_length_;

    forAll(rodsCellList_, j)
    {
      label cell_j = rodsCellList_[j];
      cell_volumes[cell_j] -= 1.0 / 4.0 * rod_section_volume_;

      if (cell_volumes[cell_j] < 0.0)
      {
        FatalErrorIn("void Foam::subChannelMesh::modify_mesh()")
                  << "Cell Volume is negative !"
                  << exit(FatalError);
      }
    }
  }


  /*--------------------------------------------------------------------------*\
  | internal faces                                                             |
  \*--------------------------------------------------------------------------*/

  surfaceVectorField * SfPtr_ = get_SfPtr();
  surfaceVectorField & faces_area_vector = *SfPtr_;

  surfaceScalarField * magSfPtr_ = get_magSfPtr();
  surfaceScalarField & faces_area_mag = *magSfPtr_;

  forAll(rodsSectionList_, i)
  {
    label edge_i = (*rodEdgesListPtr_)[i];
    labelList rodsFaceList_ = edgeFaces()[edge_i];

    scalar rod_outer_diameter_ = RodOuterDiameter_;
    scalar rod_section_length_ = get_rod_section_length(edge_i);
    scalar rod_section_area_   = 1.0 / 2.0 * rod_outer_diameter_ * rod_section_length_;

    forAll(rodsFaceList_, j)
    {
      label face_j = rodsFaceList_[j];

      if (face_j < faces_area_vector.size())
      {
        if(check_face_vertical_to_z(faces_area_vector[face_j]))
        {
          vector face_area_vector = faces_area_vector[face_j];

          faces_area_vector[face_j] = vector_remove_a_scalar(face_area_vector, rod_section_area_);

          faces_area_mag[face_j] -= rod_section_area_;

          if (faces_area_mag[face_j] < 0.0)
          {
            FatalErrorIn("void Foam::subChannelMesh::modify_mesh()")
                      << "Face Area is negative !"
                      << exit(FatalError);
          }
        }
      }
      else if (Pstream::parRun())
      {
        label find_face_times = 0;

        forAll(this->boundaryMesh(), boundaryMeshI)
        {
          label startFace_ = this->boundaryMesh()[boundaryMeshI].start();
          label nFaces_ = this->boundaryMesh()[boundaryMeshI].size();

          word boundary_group = this->boundaryMesh()[boundaryMeshI].inGroups()[0];

          if ( (face_j >= startFace_) && (face_j < (startFace_+nFaces_)) && (boundary_group == "processor") )
          {
            find_face_times += 1;

            label boundary_i = boundaryMeshI;
            label boundary_face_i = face_j - startFace_;

            if(check_face_vertical_to_z(faces_area_vector.boundaryField()[boundary_i][boundary_face_i]))
            {
              vector face_area_vector = faces_area_vector.boundaryField()[boundary_i][boundary_face_i];

              faces_area_vector.boundaryField()[boundary_i][boundary_face_i] =
                vector_remove_a_scalar(face_area_vector, rod_section_area_);

              faces_area_mag.boundaryField()[boundary_i][boundary_face_i] -= rod_section_area_;

              if (faces_area_mag.boundaryField()[boundary_i][boundary_face_i] < 0.0)
              {
                FatalErrorIn("void Foam::subChannelMesh::modify_mesh()")
                          << "Face Area is negative !"
                          << exit(FatalError);
              }
            }
          }
        }

        if (find_face_times != 1)
        {
          FatalErrorIn("void Foam::subChannelMesh::modify_mesh()")
                    << "Boundary face in processor group is not found just once !"
                    << exit(FatalError);
        }
      }
      else
      {
        FatalErrorIn("void Foam::subChannelMesh::modify_mesh()")
                  << "Face is on the boundary and it is serial run !"
                  << exit(FatalError);
      }
    }
  }


  scalar pi_ = 3.14159;
  scalar rod_outer_diameter_ = RodOuterDiameter_;
  scalar rod_cross_section_area_ = 1.0 / 4.0 * pi_ * ::sqr(rod_outer_diameter_);

  labelList pointListCheck(this->nPoints(), 1);

  forAll(rodsSectionList_, i)
  {
    label edge_i = (*rodEdgesListPtr_)[i];

    label startLabel = this->edges()[edge_i].start();
    label endLabel   = this->edges()[edge_i].end();

    labelList pointFacesList = this->pointFaces()[startLabel];

    if (pointListCheck[startLabel])
    {
      forAll(pointFacesList, j)
      {
        label face_j = pointFacesList[j];

        if (face_j < faces_area_vector.size())
        {
          if (!check_face_vertical_to_z(faces_area_vector[face_j]))
          {

            vector face_area_vector = faces_area_vector[face_j];

            faces_area_vector[face_j] =
              vector_remove_a_scalar(face_area_vector, 1.0 / 4.0 * rod_cross_section_area_);

            faces_area_mag[face_j] -= 1.0 / 4.0 * rod_cross_section_area_;

            if (faces_area_mag[face_j] < 0.0)
            {
              FatalErrorIn("void Foam::subChannelMesh::modify_mesh()")
                        << "Face Area is negative !"
                        << exit(FatalError);
            }
          }
        }
        else if (Pstream::parRun())
        {
          label find_face_times = 0;

          forAll(this->boundaryMesh(), boundaryMeshI)
          {
            label startFace_ = this->boundaryMesh()[boundaryMeshI].start();
            label nFaces_ = this->boundaryMesh()[boundaryMeshI].size();

            word boundary_group = this->boundaryMesh()[boundaryMeshI].inGroups()[0];

            if ( (face_j >= startFace_) && (face_j < (startFace_+nFaces_)) && (boundary_group == "processor") )
            {
              find_face_times += 1;

              label boundary_i = boundaryMeshI;
              label boundary_face_i = face_j - startFace_;

              if(!check_face_vertical_to_z(faces_area_vector.boundaryField()[boundary_i][boundary_face_i]))
              {
                vector face_area_vector = faces_area_vector.boundaryField()[boundary_i][boundary_face_i];

                faces_area_vector.boundaryField()[boundary_i][boundary_face_i] =
                  vector_remove_a_scalar(face_area_vector, 1.0 / 4.0 *rod_cross_section_area_);

                faces_area_mag.boundaryField()[boundary_i][boundary_face_i] -= 1.0 / 4.0 * rod_cross_section_area_;

                if (faces_area_mag.boundaryField()[boundary_i][boundary_face_i] < 0.0)
                {
                  FatalErrorIn("void Foam::subChannelMesh::modify_mesh()")
                            << "Face Area is negative !"
                            << exit(FatalError);
                }
              }
            }
          }

          if (find_face_times > 1)
          {
            FatalErrorIn("void Foam::subChannelMesh::modify_mesh()")
                      << "Boundary face in processor group is not found just once !"
                      << exit(FatalError);
          }
        }
      }
      pointListCheck[startLabel] = 0;
    }


    pointFacesList = this->pointFaces()[endLabel];

    if (pointListCheck[endLabel])
    {
      forAll(pointFacesList, j)
      {
        label face_j = pointFacesList[j];

        if (face_j < faces_area_vector.size())
        {
          if (!check_face_vertical_to_z(faces_area_vector[face_j]))
          {

            vector face_area_vector = faces_area_vector[face_j];

            faces_area_vector[face_j] =
              vector_remove_a_scalar(face_area_vector, 1.0 / 4.0 * rod_cross_section_area_);

            faces_area_mag[face_j] -= 1.0 / 4.0 * rod_cross_section_area_;

            if (faces_area_mag[face_j] < 0.0)
            {
              FatalErrorIn("void Foam::subChannelMesh::modify_mesh()")
                        << "Face Area is negative !"
                        << exit(FatalError);
            }
          }
        }
        else if (Pstream::parRun())
        {
          label find_face_times = 0;

          forAll(this->boundaryMesh(), boundaryMeshI)
          {
            label startFace_ = this->boundaryMesh()[boundaryMeshI].start();
            label nFaces_ = this->boundaryMesh()[boundaryMeshI].size();

            word boundary_group = this->boundaryMesh()[boundaryMeshI].inGroups()[0];

            if ( (face_j >= startFace_) && (face_j < (startFace_+nFaces_)) && (boundary_group == "processor") )
            {
              find_face_times += 1;

              label boundary_i = boundaryMeshI;
              label boundary_face_i = face_j - startFace_;

              if(!check_face_vertical_to_z(faces_area_vector.boundaryField()[boundary_i][boundary_face_i]))
              {
                vector face_area_vector = faces_area_vector.boundaryField()[boundary_i][boundary_face_i];

                faces_area_vector.boundaryField()[boundary_i][boundary_face_i] =
                  vector_remove_a_scalar(face_area_vector, 1.0 / 4.0 *rod_cross_section_area_);

                faces_area_mag.boundaryField()[boundary_i][boundary_face_i] -= 1.0 / 4.0 * rod_cross_section_area_;

                if (faces_area_mag.boundaryField()[boundary_i][boundary_face_i] < 0.0)
                {
                  FatalErrorIn("void Foam::subChannelMesh::modify_mesh()")
                            << "Face Area is negative !"
                            << exit(FatalError);
                }
              }
            }
          }

          if (find_face_times > 1)
          {
            FatalErrorIn("void Foam::subChannelMesh::modify_mesh()")
                      << "Boundary face in processor group is not found just once !"
                      << exit(FatalError);
          }
        }
      }
      pointListCheck[endLabel] = 0;
    }

  }


  /*--------------------------------------------------------------------------*\
  | boundary                                                                   |
  \*--------------------------------------------------------------------------*/

  word boundaryName_("INLET");
  label boundary_i = get_boundary_label_by_name(boundaryName_);

  label boundaryStart_ = 0;
  label boundarySize_  = 0;

  get_boundary_information_by_name(boundaryName_, boundaryStart_, boundarySize_);

  pointListCheck.clear();
  pointListCheck.setSize(this->nPoints(), 1);


  forAll(rodsSectionList_, i)
  {

    label edge_i = (*rodEdgesListPtr_)[i];

    label startLabel = this->edges()[edge_i].start();
    label endLabel   = this->edges()[edge_i].end();

    labelList pointFacesList = this->pointFaces()[startLabel];

    if (pointListCheck[startLabel])
    {
      forAll(pointFacesList, faceI)
      {
        if ((pointFacesList[faceI] >= boundaryStart_) && (pointFacesList[faceI] < (boundaryStart_ + boundarySize_)))
        {
          label boundary_face_i = pointFacesList[faceI] - boundaryStart_;

          vector face_area_vector = faces_area_vector.boundaryField()[boundary_i][boundary_face_i];

          faces_area_vector.boundaryField()[boundary_i][boundary_face_i] =
            vector_remove_a_scalar(face_area_vector, 1.0 / 4.0 * rod_cross_section_area_);

          faces_area_mag.boundaryField()[boundary_i][boundary_face_i] -= 1.0 / 4.0 * rod_cross_section_area_;

          if (faces_area_mag.boundaryField()[boundary_i][boundary_face_i] < 0.0)
          {
            FatalErrorIn("void Foam::subChannelMesh::modify_mesh()")
                      << "Face Area is negative !"
                      << exit(FatalError);
          }
        }
      }
      pointListCheck[startLabel] = 0;
    }


    pointFacesList = this->pointFaces()[endLabel];

    if (pointListCheck[endLabel])
    {
      forAll(pointFacesList, faceI)
      {
        if ((pointFacesList[faceI] >= boundaryStart_) && (pointFacesList[faceI] < (boundaryStart_ + boundarySize_)))
        {
          label boundary_face_i = pointFacesList[faceI] - boundaryStart_;

          vector face_area_vector = faces_area_vector.boundaryField()[boundary_i][boundary_face_i];

          faces_area_vector.boundaryField()[boundary_i][boundary_face_i] =
            vector_remove_a_scalar(face_area_vector, 1.0 / 4.0 * rod_cross_section_area_);

          faces_area_mag.boundaryField()[boundary_i][boundary_face_i] -= 1.0 / 4.0 * rod_cross_section_area_;

          if (faces_area_mag.boundaryField()[boundary_i][boundary_face_i] < 0.0)
          {
            FatalErrorIn("void Foam::subChannelMesh::modify_mesh()")
                      << "Face Area is negative !"
                      << exit(FatalError);
          }
        }
      }
      pointListCheck[endLabel] = 0;
    }
  }


  boundaryName_ = "OUTLET";
  boundary_i = get_boundary_label_by_name(boundaryName_);

  boundaryStart_ = 0;
  boundarySize_  = 0;

  get_boundary_information_by_name(boundaryName_, boundaryStart_, boundarySize_);

  pointListCheck.clear();
  pointListCheck.setSize(this->nPoints(), 1);

  forAll(rodsSectionList_, i)
  {

    label edge_i = (*rodEdgesListPtr_)[i];

    label startLabel = this->edges()[edge_i].start();
    label endLabel   = this->edges()[edge_i].end();

    labelList pointFacesList = this->pointFaces()[startLabel];

    if (pointListCheck[startLabel])
    {
      forAll(pointFacesList, faceI)
      {
        if ((pointFacesList[faceI] >= boundaryStart_) && (pointFacesList[faceI] < (boundaryStart_ + boundarySize_)))
        {
          label boundary_face_i = pointFacesList[faceI] - boundaryStart_;

          vector face_area_vector = faces_area_vector.boundaryField()[boundary_i][boundary_face_i];

          faces_area_vector.boundaryField()[boundary_i][boundary_face_i] =
            vector_remove_a_scalar(face_area_vector, 1.0 / 4.0 * rod_cross_section_area_);

          faces_area_mag.boundaryField()[boundary_i][boundary_face_i] -= 1.0 / 4.0 * rod_cross_section_area_;

          if (faces_area_mag.boundaryField()[boundary_i][boundary_face_i] < 0.0)
          {
            FatalErrorIn("void Foam::subChannelMesh::modify_mesh()")
                      << "Face Area is negative !"
                      << exit(FatalError);
          }
        }
      }
      pointListCheck[startLabel] = 0;
    }


    pointFacesList = this->pointFaces()[endLabel];

    if (pointListCheck[endLabel])
    {
      forAll(pointFacesList, faceI)
      {
        if ((pointFacesList[faceI] >= boundaryStart_) && (pointFacesList[faceI] < (boundaryStart_ + boundarySize_)))
        {
          label boundary_face_i = pointFacesList[faceI] - boundaryStart_;

          vector face_area_vector = faces_area_vector.boundaryField()[boundary_i][boundary_face_i];

          faces_area_vector.boundaryField()[boundary_i][boundary_face_i] =
            vector_remove_a_scalar(face_area_vector, 1.0 / 4.0 * rod_cross_section_area_);

          faces_area_mag.boundaryField()[boundary_i][boundary_face_i] -= 1.0 / 4.0 * rod_cross_section_area_;

          if (faces_area_mag.boundaryField()[boundary_i][boundary_face_i] < 0.0)
          {
            FatalErrorIn("void Foam::subChannelMesh::modify_mesh()")
                      << "Face Area is negative !"
                      << exit(FatalError);
          }
        }
      }
      pointListCheck[endLabel] = 0;
    }
  }

  mesh_modified_or_not = true;

  calc_flow_area();
  calc_hydraulic_diameter();
}


void
Foam::subChannelMesh::modify_mesh_for_core_edge()
{

  if (!rodsBuildOrNot)
  {
    FatalErrorIn("void Foam::subChannelMesh::modify_mesh()")
              << "rods should be built first !"
              << exit(FatalError);
  }

  labelList rodsSectionList_ = *rodEdgesListOnCoreEdgePtr_;


  surfaceVectorField * SfPtr_ = get_SfPtr();
  surfaceVectorField & faces_area_vector = *SfPtr_;

  surfaceScalarField * magSfPtr_ = get_magSfPtr();
  surfaceScalarField & faces_area_mag = *magSfPtr_;


  word boundaryName_("WALL");
  label boundary_i = get_boundary_label_by_name(boundaryName_);

  label boundaryStart_ = 0;
  label boundarySize_  = 0;

  get_boundary_information_by_name(boundaryName_, boundaryStart_, boundarySize_);


  forAll(rodsSectionList_, i)
  {
    label edge_i = (*rodEdgesListOnCoreEdgePtr_)[i];
    labelList rodsFaceList_ = edgeFaces()[edge_i];

    scalar rod_outer_diameter_ = RodOuterDiameter_;
    scalar rod_section_length_ = get_rod_section_length(edge_i);
    scalar rod_section_area_   = 1.0 / 2.0 * rod_outer_diameter_ *rod_section_length_;

    forAll(rodsFaceList_, j)
    {
      label face_j = rodsFaceList_[j];

      if (face_j < this->nInternalFaces())
      {
        vector face_area_vector = faces_area_vector[face_j];

        if (mag(face_area_vector) > rod_section_area_)
        {
          faces_area_vector[face_j] = vector_remove_a_scalar(face_area_vector, rod_section_area_);

          faces_area_mag[face_j] -= rod_section_area_;
        }

        if (faces_area_mag[face_j] < 0.0)
        {
          FatalErrorIn("void Foam::subChannelMesh::modify_mesh_for_core_edge()")
                    << "Face Area is negative !"
                    << exit(FatalError);
        }
      }
      else if ( (face_j >= boundaryStart_) && (face_j < (boundaryStart_ + boundarySize_)) )
      {
        label boundary_face_i = face_j - boundaryStart_;

        vector face_area_vector = faces_area_vector.boundaryField()[boundary_i][boundary_face_i];

        if (mag(face_area_vector) > rod_section_area_)
        {
        faces_area_vector.boundaryField()[boundary_i][boundary_face_i] =
          vector_remove_a_scalar(face_area_vector, rod_section_area_);

        faces_area_mag.boundaryField()[boundary_i][boundary_face_i] -= rod_section_area_;
        }

        if (faces_area_mag.boundaryField()[boundary_i][boundary_face_i] < 0.0)
        {
          FatalErrorIn("void Foam::subChannelMesh::modify_mesh_for_core_edge()")
                    << "Face Area is negative !"
                    << exit(FatalError);
        }
      }
    }
  }

}


void
Foam::subChannelMesh::modify_mesh_for_assembly_edge()
{
  if (!rodsBuildOrNot)
  {
    FatalErrorIn("void Foam::subChannelMesh::modify_mesh()")
              << "rods should be built first !"
              << exit(FatalError);
  }

  labelList rodsSectionList_ = *rodEdgesListOnAssemblyEdgePtr_;


  surfaceVectorField * SfPtr_ = get_SfPtr();
  surfaceVectorField & faces_area_vector = *SfPtr_;

  surfaceScalarField * magSfPtr_ = get_magSfPtr();
  surfaceScalarField & faces_area_mag = *magSfPtr_;

  forAll(rodsSectionList_, i)
  {
    label edge_i = (*rodEdgesListOnAssemblyEdgePtr_)[i];
    labelList rodsFaceList_ = edgeFaces()[edge_i];

    scalar rod_outer_diameter_ = RodOuterDiameter_;
    scalar rod_section_length_ = get_rod_section_length(edge_i);
    scalar rod_section_area_   = 1.0 / 2.0 * rod_outer_diameter_ *rod_section_length_;

    forAll(rodsFaceList_, j)
    {
      label face_j = rodsFaceList_[j];

      if (face_j < this->nInternalFaces())
      {
        vector face_area_vector = faces_area_vector[face_j];

        if (mag(face_area_vector) > rod_section_area_)
        {
          faces_area_vector[face_j] = vector_remove_a_scalar(face_area_vector, rod_section_area_);

          faces_area_mag[face_j] -= rod_section_area_;
        }

        if (faces_area_mag[face_j] < 0.0)
        {
          FatalErrorIn("void Foam::subChannelMesh::modify_mesh_for_core_edge()")
                    << "Face Area is negative !"
                    << exit(FatalError);
        }
      }
      else
      {
        FatalErrorIn("void Foam::subChannelMesh::modify_mesh_for_core_edge()")
                  << "face_j should less than nInternalFaces ! "
                  << exit(FatalError);
      }
    }
  }


}


void
Foam::subChannelMesh::createRods()
{
  if (!rodsBuildOrNot)
  {
    findRodEdgesList();
    findRodEdgesOnAssemblyEdgeList();
    findRodEdgesOnCoreEdgeList();
    findGapList();
    rodsBuildOrNot = true;

    build_cellRodsList();

    calc_index_assembly_and_rod();

    calcRodEdgesLengthList();
    calcRodEdgesCenterList();

    calc_core_z_position();

  }
  else
  {
    FatalErrorIn("void Foam::subChannelMesh::createRods()")
                << "rods shouldn't be built again !"
                << exit(FatalError);
  }
}


const labelList &
Foam::subChannelMesh::get_rodEdgesList() const
{

  if (!rodEdgesListPtr_)
  {
    FatalErrorIn("void Foam::subChannelMesh::get_rodEdgesList()")
                << "Can't find rods !"
                << exit(FatalError);
  }

  return *rodEdgesListPtr_;
}


const labelList &
Foam::subChannelMesh::get_rodEdgesListOnAssemblyEdge() const
{

  if (!rodEdgesListOnAssemblyEdgePtr_)
  {
    FatalErrorIn("void Foam::subChannelMesh::get_rodEdgesListOnAssemblyEdge()")
                << "Can't find rods !"
                << exit(FatalError);
  }

  return *rodEdgesListOnAssemblyEdgePtr_;
}


const labelList &
Foam::subChannelMesh::get_rodEdgesListOnCoreEdge() const
{

  if (!rodEdgesListOnCoreEdgePtr_)
  {
    FatalErrorIn("void Foam::subChannelMesh::get_rodEdgesListOnCoreEdge()")
                << "Can't find rods !"
                << exit(FatalError);
  }

  return *rodEdgesListOnCoreEdgePtr_;
}


const labelList &
Foam::subChannelMesh::get_gapList() const
{

  if (!gapListPtr_)
  {
    FatalErrorIn("void Foam::subChannelMesh::get_gapList()")
                << "Can't find gaps !"
                << exit(FatalError);
  }

  return *gapListPtr_;
}


const scalarList &
Foam::subChannelMesh::get_rodEdgesLengthList() const
{

  if (!rodEdgesLengthListPtr_)
  {
    FatalErrorIn("void Foam::subChannelMesh::get_rodEdgesLengthList()")
                << "Can't find rods edge length !"
                << exit(FatalError);
  }

  return *rodEdgesLengthListPtr_;
}


const pointField &
Foam::subChannelMesh::get_rodEdgesCenterList() const
{
  if (!rodEdgesCenterListPtr_)
  {
    FatalErrorIn("void Foam::subChannelMesh::get_rodEdgesCenterList()")
                << "Can't find rods edge center !"
                << exit(FatalError);
  }

  return *rodEdgesCenterListPtr_;
}


scalar
Foam::subChannelMesh::get_rod_edge_z_max(label edge_i) const
{
  scalar z_max_return = -1.0 * GREAT;

  label startLabel = this->edges()[edge_i].start();
  vector startPoint = this->points()[startLabel];
  z_max_return = (startPoint.z() > z_max_return) ? startPoint.z() : z_max_return;

  label endLabel = this->edges()[edge_i].end();
  vector endPoint = this->points()[endLabel];
  z_max_return = (endPoint.z() > z_max_return) ? endPoint.z() : z_max_return;

  return z_max_return;
}


scalar
Foam::subChannelMesh::get_rod_edge_z_min(label edge_i) const
{
  scalar z_min_return = GREAT;

  label startLabel = this->edges()[edge_i].start();
  vector startPoint = this->points()[startLabel];
  z_min_return = (startPoint.z() < z_min_return) ? startPoint.z() : z_min_return;

  label endLabel = this->edges()[edge_i].end();
  vector endPoint = this->points()[endLabel];
  z_min_return = (endPoint.z() < z_min_return) ? endPoint.z() : z_min_return;

  return z_min_return;
}


void
Foam::subChannelMesh::findRodEdgesList() const
{

  const labelListList & fcsEdges = this->faceEdges();
  labelList edgesListCheck(this->nEdges(), 1);


  forAll(this->boundaryMesh(), boundaryMeshI)
  {
    label startFace_ = this->boundaryMesh()[boundaryMeshI].start();
    label nFaces_ = this->boundaryMesh()[boundaryMeshI].size();

    if (this->boundaryMesh()[boundaryMeshI].inGroups()[0] != "processor")
    {
      for (label i = startFace_; i < startFace_ + nFaces_; i++)
      {
        label edgeI;
        forAll(fcsEdges[i], j)
        {
          edgeI = fcsEdges[i][j];
          if (edgesListCheck[edgeI])
          {
            edgesListCheck[edgeI] = 0;
          }
        }
      }
    }
  }


  forAll(edgesListCheck, edgeI)
  {
    if (edgesListCheck[edgeI])
    {
      label startLabel = this->edges()[edgeI].start();
      vector startPoint = this->points()[startLabel];

      label endLabel = this->edges()[edgeI].end();
      vector endPoint = this->points()[endLabel];

      vector edgeVector = endPoint - startPoint;

      if (mag(edgeVector.z()) < SMALL)
      {
        edgesListCheck[edgeI] = 0;
      }
    }
  }


  forAll(edgesListCheck, edgeI)
  {
    if (edgesListCheck[edgeI])
    {
      label startLabel = this->edges()[edgeI].start();
      vector startPoint = this->points()[startLabel];

      if (checkOnAssemblyEdge(startPoint))
      {
        edgesListCheck[edgeI] = 0;
      }
    }
  }


  label n = 0;
  forAll(edgesListCheck, edgeI)
  {
    if (edgesListCheck[edgeI])
    {
      n++;
    }
  }
  rodEdgesListPtr_ = new labelList(n, 0);

  n = 0;
  forAll(edgesListCheck, edgeI)
  {
    if (edgesListCheck[edgeI] == 1)
    {
      (*rodEdgesListPtr_)[n] = edgeI;
      n++;
    }
  }
}


void
Foam::subChannelMesh::findRodEdgesOnAssemblyEdgeList() const
{

  const labelListList & fcsEdges = this->faceEdges();
  labelList edgesListCheck(this->nEdges(), 1);

  forAll(this->boundaryMesh(), boundaryMeshI)
  {
    label startFace_ = this->boundaryMesh()[boundaryMeshI].start();
    label nFaces_ = this->boundaryMesh()[boundaryMeshI].size();

    if (this->boundaryMesh()[boundaryMeshI].inGroups()[0] != "processor")
    {
      for (label i = startFace_; i < startFace_ + nFaces_; i++)
      {
        label edgeI;
        forAll(fcsEdges[i], j)
        {
          edgeI = fcsEdges[i][j];
          if (edgesListCheck[edgeI])
          {
            edgesListCheck[edgeI] = 0;
          }
        }
      }
    }
  }


  forAll(edgesListCheck, edgeI)
  {
    if (edgesListCheck[edgeI])
    {
      label startLabel = this->edges()[edgeI].start();
      vector startPoint = this->points()[startLabel];

      label endLabel = this->edges()[edgeI].end();
      vector endPoint = this->points()[endLabel];

      vector edgeVector = endPoint - startPoint;

      if (mag(edgeVector.z()) < SMALL)
      {
        edgesListCheck[edgeI] = 0;
      }
    }
  }


  forAll(edgesListCheck, edgeI)
  {
    if (edgesListCheck[edgeI])
    {
      label startLabel = this->edges()[edgeI].start();
      vector startPoint = this->points()[startLabel];

      if (!checkOnAssemblyEdge(startPoint))
      {
        edgesListCheck[edgeI] = 0;
      }
    }
  }

  // Reorder the edges and build it

  label n = 0;
  forAll(edgesListCheck, edgeI)
  {
    if (edgesListCheck[edgeI])
    {
      n++;
    }
  }
  rodEdgesListOnAssemblyEdgePtr_ = new labelList(n, 0);

  n = 0;
  forAll(edgesListCheck, edgeI)
  {
    if (edgesListCheck[edgeI] == 1)
    {
      (*rodEdgesListOnAssemblyEdgePtr_)[n] = edgeI;
      n++;
    }
  }
}


void
Foam::subChannelMesh::findRodEdgesOnCoreEdgeList() const
{

  const labelListList & fcsEdges = this->faceEdges();
  labelList edgesListCheck(this->nEdges(), 1);

  forAll(this->boundaryMesh(), boundaryMeshI)
  {
    label startFace_ = this->boundaryMesh()[boundaryMeshI].start();
    label nFaces_ = this->boundaryMesh()[boundaryMeshI].size();

    if (this->boundaryMesh()[boundaryMeshI].inGroups()[0] != "processor")
    {
      for (label i = startFace_; i < startFace_ + nFaces_; i++)
      {
        label edgeI;
        forAll(fcsEdges[i], j)
        {
          edgeI = fcsEdges[i][j];
          if (edgesListCheck[edgeI])
          {
            edgesListCheck[edgeI] = 0;
          }
        }
      }
    }
  }

  forAll(edgesListCheck, i)
  {
    if(edgesListCheck[i])
    {
      edgesListCheck[i] = 0;
    }
    else
    {
      edgesListCheck[i] = 1;
    }
  }

  forAll(edgesListCheck, edgeI)
  {
    if (edgesListCheck[edgeI])
    {
      label startLabel = this->edges()[edgeI].start();
      vector startPoint = this->points()[startLabel];

      label endLabel = this->edges()[edgeI].end();
      vector endPoint = this->points()[endLabel];

      vector edgeVector = endPoint - startPoint;

      if (mag(edgeVector.z()) < SMALL)
      {
        edgesListCheck[edgeI] = 0;
      }
    }
  }

  // Reorder the edges and build

  label n = 0;
  forAll(edgesListCheck, edgeI)
  {
    if (edgesListCheck[edgeI])
    {
      n++;
    }
  }
  rodEdgesListOnCoreEdgePtr_ = new labelList(n, 0);

  n = 0;
  forAll(edgesListCheck, edgeI)
  {
    if (edgesListCheck[edgeI] == 1)
    {
      (*rodEdgesListOnCoreEdgePtr_)[n] = edgeI;
      n++;
    }
  }
}


void
Foam::subChannelMesh::findGapList()
{

  surfaceVectorField * SfPtr_ = get_SfPtr();
  surfaceVectorField & faces_area_vector = *SfPtr_;

  label gap_number_ = 0;

  forAll(faces_area_vector, faceI)
  {
    if (check_face_vertical_to_z(faces_area_vector[faceI]))
    {
      gap_number_ += 1;
    }
  }


  gapListPtr_ = new labelList(gap_number_, 0);

  label n = 0;
  forAll(faces_area_vector, faceI)
  {
    if (check_face_vertical_to_z(faces_area_vector[faceI]))
    {
      (*gapListPtr_)[n] = faceI;
      n++;
    }
  }

  if (!(n == gap_number_))
  {
    FatalErrorIn("void Foam::subChannelMesh::findGapList()")
              << "n should equal to gap_number_ !"
              << exit(FatalError);
  }
}


void
Foam::subChannelMesh::build_cellRodsList()
{
  if (!rodEdgesListPtr_)
  {
    FatalErrorIn("void Foam::subChannelMesh::build_cellRodsList()")
              << "rodEdgesListPtr_ should be built first !"
              << exit(FatalError);
  }

  if (cellRodsNumberList_)
  {
    FatalErrorIn("void Foam::subChannelMesh::build_cellRodsList()")
              << "cellRodsNumberList_ shouldn't be built again !"
              << exit(FatalError);
  }

  if (cellRodsList_)
  {
    FatalErrorIn("void Foam::subChannelMesh::build_cellRodsList()")
              << "cellRodsList_ shouldn't be built again !"
              << exit(FatalError);
  }


  volVectorField * CPtr_ = get_CPtr();
  volVectorField & cell_center_vector = *CPtr_;

  cellRodsNumberList_ = new labelList(cell_center_vector.size(), 0);

  cellRodsList_ =
    new labelListList(cell_center_vector.size(), labelList(4, -1));

  forAll((*rodEdgesListPtr_), i)
  {
    label edge_i = (*rodEdgesListPtr_)[i];

    labelList rodsCellList_ = edgeCells()[edge_i];

    forAll(rodsCellList_, j)
    {
      label cell_j = rodsCellList_[j];

      (*cellRodsNumberList_)[cell_j] += 1;

      forAll((*cellRodsList_)[cell_j], k)
      {
        if ((*cellRodsList_)[cell_j][k] < 0)
        {
          (*cellRodsList_)[cell_j][k] = edge_i;
          break;
        }
      }
    }
  }

  forAll((*cellRodsNumberList_), i)
  {
    if ((*cellRodsNumberList_)[i] > 4)
    {
      FatalErrorIn("void Foam::subChannelMesh::build_cellRodsList()")
                << "cellRodsNumberList_[i] is more than 4 !"
                << exit(FatalError);
    }
  }
}


void
Foam::subChannelMesh::calc_index_assembly_and_rod()
{
  if (rod_assembly_index_X_Y_)
  {
    FatalErrorIn("void Foam::subChannelMesh::calc_index_assembly_and_rod()")
              << "rod_assembly_index_X_Y_ can't be made again !"
              << exit(FatalError);
  }

  if (rod_rod_index_X_Y_)
  {
    FatalErrorIn("void Foam::subChannelMesh::calc_index_assembly_and_rod()")
              << "rod_rod_index_X_Y_ can't be made again !"
              << exit(FatalError);
  }

  if (!rodEdgesListPtr_)
  {
    FatalErrorIn("void Foam::subChannelMesh::calcRodEdgesLengthList()")
              << "rodEdgesListPtr_ should be built first !"
              << exit(FatalError);
  }

  rod_assembly_index_X_Y_ =
    new labelListList((*rodEdgesListPtr_).size(), labelList(2, -1));

  rod_rod_index_X_Y_ =
    new labelListList((*rodEdgesListPtr_).size(), labelList(2, -1));


  forAll((*rodEdgesListPtr_), i)
  {
    label edge_i = (*rodEdgesListPtr_)[i];
    vector edge_center = get_rod_section_center(edge_i);

    label x_found_number = 0;
    label y_found_number = 0;

    for (label x_i=1; x_i<CoreMapDimension_+1; x_i++)
    {
      if (edge_center.x() > (x_i - CenterAssemblyXIndex_ - 0.5)*AssemblyWidth_ &&
          edge_center.x() < (x_i - CenterAssemblyXIndex_ + 0.5)*AssemblyWidth_)
      {
        (*rod_assembly_index_X_Y_)[i][0] = x_i;
        x_found_number += 1;
      }
    }

    for (label y_j=1; y_j<CoreMapDimension_+1; y_j++)
    {
      if (edge_center.y() > (y_j - CenterAssemblyYIndex_ - 0.5)*AssemblyWidth_ &&
          edge_center.y() < (y_j - CenterAssemblyYIndex_ + 0.5)*AssemblyWidth_)
      {
        (*rod_assembly_index_X_Y_)[i][1] = y_j;
        y_found_number += 1;
      }
    }

    if (!((x_found_number == 1) && (y_found_number == 1)))
    {
      FatalErrorIn("void Foam::subChannelMesh::calc_index_assembly_and_rod()")
                << "Error occurs when finding assembly ! "
                << exit(FatalError);
    }
  }


  if (correct_rod_arrangement_)
  {
    forAll((*rodEdgesListPtr_), i)
    {
      label edge_i = (*rodEdgesListPtr_)[i];
      vector edge_center = get_rod_section_center(edge_i);

      label x_found_number = 0;
      label y_found_number = 0;

      label assembly_index_X = (*rod_assembly_index_X_Y_)[i][0];
      label assembly_index_Y = (*rod_assembly_index_X_Y_)[i][1];

      scalar assembly_x = (assembly_index_X - CenterAssemblyXIndex_)*AssemblyWidth_;
      scalar assembly_y = (assembly_index_Y - CenterAssemblyYIndex_)*AssemblyWidth_;

      scalar criteria = (AssemblyWidth_ - pitch_*(AssemblyRodsDimension_-1)) / 2.0;

      for (label x_i=1; x_i<AssemblyRodsDimension_+1; x_i++)
      {
        if (edge_center.x() > assembly_x + (x_i - center_rod_X_index_)*pitch_ - criteria &&
            edge_center.x() < assembly_x + (x_i - center_rod_X_index_)*pitch_ + criteria)
        {
          (*rod_rod_index_X_Y_)[i][0] = x_i;
          x_found_number += 1;
        }
      }

      for (label y_j=1; y_j<AssemblyRodsDimension_+1; y_j++)
      {
        if (edge_center.y() > assembly_y + (y_j - center_rod_Y_index_)*pitch_ - criteria &&
            edge_center.y() < assembly_y + (y_j - center_rod_Y_index_)*pitch_ + criteria)
        {
          (*rod_rod_index_X_Y_)[i][1] = y_j;
          y_found_number += 1;
        }
      }
    }
  }
  else
  {
    forAll((*rodEdgesListPtr_), i)
    {
      label edge_i = (*rodEdgesListPtr_)[i];
      vector edge_center = get_rod_section_center(edge_i);

      label x_found_number = 0;
      label y_found_number = 0;

      label assembly_index_X = (*rod_assembly_index_X_Y_)[i][0];
      label assembly_index_Y = (*rod_assembly_index_X_Y_)[i][1];

      scalar assembly_x = (assembly_index_X - CenterAssemblyXIndex_)*AssemblyWidth_;
      scalar assembly_y = (assembly_index_Y - CenterAssemblyYIndex_)*AssemblyWidth_;

      scalar rod_cell_width = AssemblyWidth_ / (AssemblyRodsDimension_ + 1);
      scalar criteria = rod_cell_width / 2.0;

      for (label x_i=1; x_i<AssemblyRodsDimension_+1; x_i++)
      {
        if (edge_center.x() > assembly_x + (x_i - center_rod_X_index_)*rod_cell_width - criteria &&
            edge_center.x() < assembly_x + (x_i - center_rod_X_index_)*rod_cell_width + criteria)
        {
          (*rod_rod_index_X_Y_)[i][0] = x_i;
          x_found_number += 1;
        }
      }

      for (label y_j=1; y_j<AssemblyRodsDimension_+1; y_j++)
      {
        if (edge_center.y() > assembly_y + (y_j - center_rod_Y_index_)*rod_cell_width - criteria &&
            edge_center.y() < assembly_y + (y_j - center_rod_Y_index_)*rod_cell_width + criteria)
        {
          (*rod_rod_index_X_Y_)[i][1] = y_j;
          y_found_number += 1;
        }
      }

      if (!((x_found_number == 1) && (y_found_number == 1)))
      {
        FatalErrorIn("void Foam::subChannelMesh::calc_index_assembly_and_rod()")
                  << "Error occurs when finding rods in assembly ! "
                  << exit(FatalError);
      }
    }
  }

  // Check the list
  forAll((*rod_assembly_index_X_Y_), i)
  {
    forAll((*rod_assembly_index_X_Y_)[i], j)
    {
      if ((*rod_assembly_index_X_Y_)[i][j] < 0)
      {
        FatalErrorIn("void Foam::subChannelMesh::calc_index_assembly_and_rod()")
                  << "*rod_assembly_index_X_Y_ has -1"
                  << exit(FatalError);
      }
    }
  }

  forAll((*rod_rod_index_X_Y_), i)
  {
    forAll((*rod_rod_index_X_Y_)[i], j)
    {
      if ((*rod_rod_index_X_Y_)[i][j] < 0)
      {
        FatalErrorIn("void Foam::subChannelMesh::calc_index_assembly_and_rod()")
                  << "*rod_rod_index_X_Y_ has -1"
                  << exit(FatalError);
      }
    }
  }
}


void
Foam::subChannelMesh::calcRodEdgesLengthList()
{
  if (rodEdgesLengthListPtr_)
  {
    FatalErrorIn("void Foam::subChannelMesh::calcRodEdgesLengthList()")
              << "rodEdgesLengthListPtr_ can't be made again !"
              << exit(FatalError);
  }

  if (!rodEdgesListPtr_)
  {
    FatalErrorIn("void Foam::subChannelMesh::calcRodEdgesLengthList()")
              << "rodEdgesListPtr_ should be built first !"
              << exit(FatalError);
  }

  rodEdgesLengthListPtr_ = new scalarList((*rodEdgesListPtr_).size(), 0.0);

  forAll((*rodEdgesListPtr_), i)
  {
    label edge_i = (*rodEdgesListPtr_)[i];

    (*rodEdgesLengthListPtr_)[i] = get_rod_section_length(edge_i);
  }
}


void
Foam::subChannelMesh::calcRodEdgesCenterList()
{
  if (rodEdgesCenterListPtr_)
  {
    FatalErrorIn("void Foam::subChannelMesh::calcRodEdgesCenterList()")
              << "rodEdgesCenterListPtr_ can't be made again !"
              << exit(FatalError);
  }

  if (!rodEdgesListPtr_)
  {
    FatalErrorIn("void Foam::subChannelMesh::calcRodEdgesCenterList()")
              << "rodEdgesListPtr_ should be built first !"
              << exit(FatalError);
  }

  rodEdgesCenterListPtr_ = new pointField((*rodEdgesListPtr_).size(), vector(0.0, 0.0, 0.0));

  forAll((*rodEdgesListPtr_), i)
  {
    label edge_i = (*rodEdgesListPtr_)[i];

    (*rodEdgesCenterListPtr_)[i] = get_rod_section_center(edge_i);
  }
}


void
Foam::subChannelMesh::calc_core_z_position()
{
  if (core_type_)
  {
    FatalErrorIn("void Foam::subChannelMesh::calc_core_z_position()")
              << "(*core_type_) shouldn't be built again !"
              << exit(FatalError);
  }

  volVectorField * CPtr_ = get_CPtr();
  volVectorField & cell_center_vector = *CPtr_;

  core_type_ = new labelList(cell_center_vector.size(), -1);

  forAll(cell_center_vector, i)
  {
    label found_z_times = 0;

    forAll(core_type_z_direction_, j)
    {
      if (j == 0 && cell_center_vector[i].z() < core_height_z_direction_[j])
      {
        (*core_type_)[i] = core_type_z_direction_[j];
        found_z_times += 1;
      }
      else if ( cell_center_vector[i].z() < core_height_z_direction_[j] &&
                cell_center_vector[i].z() > core_height_z_direction_[j-1] )
      {
        (*core_type_)[i] = core_type_z_direction_[j];
        found_z_times += 1;
      }
    }

    if (!(found_z_times == 1))
    {
      FatalErrorIn("void Foam::subChannelMesh::calc_core_z_position()")
                << " z position is not found only once !"
                << exit(FatalError);
    }
  }

  forAll((*core_type_), i)
  {
    if ((*core_type_)[i] < 0)
    {
      FatalErrorIn("void Foam::subChannelMesh::calc_core_z_position()")
                << "(*core_type_) has -1 !"
                << exit(FatalError);
    }
  }
}


bool
Foam::subChannelMesh::checkOnAssemblyEdge(vector p) const
{
  scalar criteria_ =
    AssemblyWidth_ / (AssemblyRodsDimension_ + 1) / 2.0;

  for (label i=1; i<CoreMapDimension_; i++)
  {
    if (p.x() > (i-CenterAssemblyXIndex_+0.5)*AssemblyWidth_ - criteria_ &&
        p.x() < (i-CenterAssemblyXIndex_+0.5)*AssemblyWidth_ + criteria_)
    {
      return true;
    }
  }

  for (label i=1; i<CoreMapDimension_; i++)
  {
    if (p.y() > (i-CenterAssemblyYIndex_+0.5)*AssemblyWidth_ - criteria_ &&
        p.y() < (i-CenterAssemblyYIndex_+0.5)*AssemblyWidth_ + criteria_)
    {
      return true;
    }
  }

  return false;
}


bool
Foam::subChannelMesh::check_face_vertical_to_z(vector p) const
{
  vector z_(0.0, 0.0, 1.0);

  if ( mag(p & z_) < SMALL )
  {
    return true;
  }
  else
  {
    return false;
  }
}


scalar
Foam::subChannelMesh::get_rod_section_length(label edgeI)
{
  scalar rod_section_length_ = 0.0;

  label startLabel = this->edges()[edgeI].start();
  vector startPoint = this->points()[startLabel];

  label endLabel = this->edges()[edgeI].end();
  vector endPoint = this->points()[endLabel];

  vector edgeVector = endPoint - startPoint;

  rod_section_length_ = mag(edgeVector);

  return rod_section_length_;
}


vector
Foam::subChannelMesh::get_rod_section_center(label edgeI)
{
  vector rod_section_center_(0.0, 0.0, 0.0);

  label startLabel = this->edges()[edgeI].start();
  vector startPoint = this->points()[startLabel];

  label endLabel = this->edges()[edgeI].end();
  vector endPoint = this->points()[endLabel];

  rod_section_center_.x() = (startPoint.x() + endPoint.x()) / 2.0;
  rod_section_center_.y() = (startPoint.y() + endPoint.y()) / 2.0;
  rod_section_center_.z() = (startPoint.z() + endPoint.z()) / 2.0;

  return rod_section_center_;
}


vector
Foam::subChannelMesh::vector_remove_a_scalar
(
  vector p,
  scalar s
)
{
  vector p_return(0.0, 0.0, 0.0);

  if (mag(p) < s)
  {
    FatalErrorIn("void Foam::subChannelMesh::vector_remove_a_scalar()")
              << "Magnitude of the vector is less than the scalar !"
              << exit(FatalError);
  }

  p_return.x() = p.x()*(1.0 - s / mag(p));
  p_return.y() = p.y()*(1.0 - s / mag(p));
  p_return.z() = p.z()*(1.0 - s / mag(p));

  return p_return;
}


void
Foam::subChannelMesh::get_boundary_information_by_name
(
  word  boundaryName,
  label & boundaryStart,
  label & boundarySize
)
{
  bool find_boundary = false;
  label find_times = 0;

  forAll(this->boundaryMesh(), boundaryI)
  {
    if (this->boundaryMesh()[boundaryI].name() == boundaryName)
    {
      boundaryStart = this->boundaryMesh()[boundaryI].start();
      boundarySize  = this->boundaryMesh()[boundaryI].size();

      find_boundary = true;
      find_times += 1;
    }
  }

  if (!find_boundary)
  {
    FatalErrorIn("void Foam::subChannelMesh::get_boundary_information_by_name()")
              << "Can't find boundary mesh !"
              << exit(FatalError);
  }

  if (find_times > 1)
  {
    FatalErrorIn("void Foam::subChannelMesh::get_boundary_information_by_name()")
              << "Find boundary mesh more than once !"
              << exit(FatalError);
  }
}


label
Foam::subChannelMesh::get_boundary_label_by_name
(
  word boundaryName
)
{
  label boundaryI_return = -1;

  forAll(this->boundaryMesh(), boundaryI)
  {
    if (this->boundaryMesh()[boundaryI].name() == boundaryName)
    {
      boundaryI_return = boundaryI;
    }
  }

  if (boundaryI_return < 0)
  {
    FatalErrorIn("void Foam::subChannelMesh::get_boundary_label_by_name()")
              << "Can't find boundary mesh !"
              << exit(FatalError);
  }

  return boundaryI_return;
}


void
Foam::subChannelMesh::calc_flow_area()
{
  scalar area_ = 0.0;
  label found_face_number = 0;

  surfaceVectorField * SfPtr_ = get_SfPtr();
  surfaceVectorField & faces_area_vector = *SfPtr_;

  surfaceScalarField * magSfPtr_ = get_magSfPtr();
  surfaceScalarField & faces_area_mag = *magSfPtr_;

  forAll(this->cells(), i)
  {
    label cell_i = i;
    labelList cell_faces_list = this->cells()[i];

    area_ = 0.0;
    found_face_number = 0;

    forAll(cell_faces_list, j)
    {
      label face_j = cell_faces_list[j];
      if (face_j < this->nInternalFaces())
      {
        if(!check_face_vertical_to_z(faces_area_vector[face_j]))
        {
          area_ += faces_area_mag[face_j];
          found_face_number += 1;
        }
      }
    }

    if (!(found_face_number > 0))
    {
      FatalErrorIn("void Foam::subChannelMesh::calc_flow_area()")
                << "Can't find faces !"
                << exit(FatalError);
    }
    else
    {
      flow_area_[cell_i] = area_ / found_face_number;
    }
  }
}


void
Foam::subChannelMesh::calc_hydraulic_diameter()
{

  forAll(this->cells(), i)
  {
    label cell_i = i;
    label cell_has_rod_number = 0;

    labelList cell_edges_list = this->cellEdges()[i];

    cell_has_rod_number = (*cellRodsNumberList_)[cell_i];

    if (!(cell_has_rod_number > 0))
    {
      FatalErrorIn("void Foam::subChannelMesh::calc_hydraulic_diameter()")
                << "Can't find rods !"
                << exit(FatalError);
    }
    else
    {
      scalar pi_ = 3.1415926;
      scalar perimeter_ = 1.0 / 4.0 * pi_ * RodOuterDiameter_ * cell_has_rod_number;

      if (test_case_with_assembly_box_)
      {
        if (cell_has_rod_number == 1)
        {
          perimeter_ += 2.0 * pi_ * assembly_box_corner_radius_ / 4.0 +
                        2.0 * ((AssemblyWidth_ - (AssemblyRodsDimension_ - 1) * pitch_)/2.0 - assembly_box_corner_radius_);

          flow_area_[cell_i] -= ::pow(assembly_box_corner_radius_, 2.0) -
                                pi_ * ::pow(assembly_box_corner_radius_, 2.0) / 4.0;
          if (flow_area_[cell_i] < 0.0)
          {
            FatalErrorIn("void Foam::subChannelMesh::calc_hydraulic_diameter()")
                      << "Assembly box corner subChannel flow area is less than 0.0 !"
                      << exit(FatalError);
          }
        }
        else if (cell_has_rod_number == 2)
        {
          perimeter_ += pitch_;
        }
        else if (cell_has_rod_number == 4)
        {
        }
        else
        {
          FatalErrorIn("void Foam::subChannelMesh::calc_hydraulic_diameter()")
                    << "subChannel cell has wrong rods number !"
                    << exit(FatalError);
        }
      }

      hydraulic_diameter_[cell_i] = 4.0 * flow_area_[cell_i] / perimeter_;
    }
  }
}


void
Foam::subChannelMesh::readReactorCore()
{

    IOdictionary ReactorCore
    (
     IOobject
     (
      "ReactorCore",
      this->time().constant(),
      this->time(),
      IOobject::MUST_READ,
      IOobject::NO_WRITE
     )
    );

    CoreMapDimension_ =
      readLabel(ReactorCore.lookup("CoreMapDimension"));

    AssemblyWidth_ =
      readScalar(ReactorCore.lookup("AssemblyWidth"));

    AssemblyRodsDimension_ =
      readLabel(ReactorCore.lookup("AssemblyRodsDimension"));

    CenterAssemblyXIndex_ =
      readLabel(ReactorCore.lookup("CenterAssemblyXIndex"));

    CenterAssemblyYIndex_ =
      readLabel(ReactorCore.lookup("CenterAssemblyYIndex"));

    center_rod_X_index_ = (AssemblyRodsDimension_) / 2.0;

    center_rod_Y_index_ = (AssemblyRodsDimension_) / 2.0;


    RodOuterDiameter_ =
      readScalar(ReactorCore.lookup("RodOuterDiameter"));

    pitch_ =
      readScalar(ReactorCore.lookup("pitch"));

    correct_rod_arrangement_ =
      ReactorCore.lookupOrDefault<Switch>("CorrectRodArrangement", Switch("no"));

    test_case_with_assembly_box_ =
      ReactorCore.lookupOrDefault<Switch>("TestCaseWithAssemblyBox", Switch("no"));

    if (test_case_with_assembly_box_)
    {
      assembly_box_corner_radius_ =
        ReactorCore.lookupOrDefault<scalar>("AssemblyBoxCornerRadius", scalar(0.0));

      if (assembly_box_corner_radius_ > (AssemblyWidth_ - (AssemblyRodsDimension_ - 1)*pitch_) / 2.0)
      {
        FatalErrorIn("void Foam::subChannelMesh::readReactorCore()")
                  << "assembly_box_corner_radius_ should less than rod outer radius !"
                  << exit(FatalError);
      }
    }

    mixing_coefficient_ =
      readScalar(ReactorCore.lookup("MixingCoefficient"));

    core_assembly_distribution_ =
      labelListList(ReactorCore.lookup("CoreAssemblyDistribution"));
    Foam::reverse(core_assembly_distribution_);

    assembly_rod_layout_ =
      labelListList(ReactorCore.lookup("AssemblyRodLayout"));
    Foam::reverse(assembly_rod_layout_);

    core_assembly_power_distribution_ =
      scalarListList(ReactorCore.lookup("CoreAssemblyPowerDistributionFactor"));
    Foam::reverse(core_assembly_power_distribution_);

    rod_pin_power_distribution_factor_ =
      scalarListList(ReactorCore.lookup("RodPinPowerDistributionFactor"));
    Foam::reverse(rod_pin_power_distribution_factor_);

    core_type_z_direction_ =
      labelList(ReactorCore.lookup("Core_type_z_direction"));

    core_height_z_direction_ =
      scalarList(ReactorCore.lookup("Core_height_z_direction"));
}


void
Foam::subChannelMesh::buildSubChannelMesh()
{
  if(!subChannelMeshBuildOrNot)
  {
    // Build fvMesh Geometry
    this->V();
    this->Sf();
    this->magSf();
    this->C();
    this->Cf();

    subChannelMeshBuildOrNot = true;
  }
  else
  {
    FatalErrorIn("void Foam::subChannelMesh::buildSubChannelMesh()")
              << "subChannelMesh can't be built again !"
              << exit(FatalError);
  }
}


void
Foam::subChannelMesh::test_function()
{
  // This is a test function...
  Info << "This is a test function !" << nl << endl;

  volScalarField::DimensionedInternalField * vptr = get_VPtr();
  volScalarField::DimensionedInternalField & cell_volumes = *vptr;

  scalar face1 = 0.0;
  scalar face2 = 0.0;

  forAll(Sf(), i)
  {
    face1 += mag(Sf()[i]);
  }

  forAll(magSf(), i)
  {
    face2 += magSf()[i];
  }

  Info << "scalar1 = " << face1 << endl;
  Info << "scalar2 = " << face2 << endl;


  word si("INLET");
  label si1 = 0;
  label si2 = 0;
  Info << "si1 = " << si1 << endl;
  Info << "si2 = " << si2 << endl;
  get_boundary_information_by_name(si, si1, si2);
  Info << "si1 = " << si1 << endl;
  Info << "si2 = " << si2 << endl;


  word g1("INLET");
  Info << get_boundary_label_by_name(g1) << endl;

  g1 = "OUTLET";
  Info << get_boundary_label_by_name(g1) << endl;

  g1 = "WALL";
  Info << get_boundary_label_by_name(g1) << endl;

  Info << "~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
  Info << hydraulic_diameter_.size() << endl;
  Info << flow_area_.size() << endl;
  Info << "~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
}


void
Foam::subChannelMesh::output_mesh_information()
{

  Info << "Mesh Information: " << nl << endl;

  volScalarField::DimensionedInternalField * vptr = get_VPtr();
  if (!vptr)
  {
    FatalErrorIn("void Foam::subChannelMesh::output_mesh_information()")
                << "Error message ...."
                << exit(FatalError);
  }

  volScalarField::DimensionedInternalField & cell_volumes = *vptr;
  scalar totalCellVolume_ = 0.0;
  label  cellVolumeNumber_ = cell_volumes.size();

  forAll(cell_volumes, i)
  {
    totalCellVolume_ += cell_volumes[i];
  }

  Info << "  Cell Volume informations :" << endl;
  Info << "  Cell Volume Number : " << cellVolumeNumber_ << endl;
  Info << "  Total Cell Volume  : " << totalCellVolume_ << "  (m^3)" << nl << endl;


  // Output

  surfaceVectorField * SfPtr_ = get_SfPtr();
  surfaceVectorField & faces_area_vector = *SfPtr_;

  surfaceScalarField * magSfPtr_ = get_magSfPtr();
  surfaceScalarField & faces_area_mag = *magSfPtr_;

  label face_parallel_to_rod_number_ = 0;
  scalar total_face_area_parallel_to_rod_ = 0.0;

  forAll(faces_area_vector, faceI)
  {
    if(check_face_vertical_to_z(faces_area_vector[faceI]))
    {
      face_parallel_to_rod_number_ += 1;
      total_face_area_parallel_to_rod_ += faces_area_mag[faceI];
    }
  }

  Info << "  Faces parallel to rod informations :" << endl;
  Info << "  Faces parallel to rod Number : " << face_parallel_to_rod_number_ << endl;
  Info << "  Total face area : " << total_face_area_parallel_to_rod_ << "  (m^2)" << nl << endl;


  label face_vertical_to_rod_number_ = 0;
  scalar total_face_area_vertical_to_rod_ = 0.0;

  forAll(faces_area_vector, faceI)
  {
    if(!check_face_vertical_to_z(faces_area_vector[faceI]))
    {
      face_vertical_to_rod_number_ += 1;
      total_face_area_vertical_to_rod_ += faces_area_mag[faceI];
    }
  }

  Info << "  Faces vertical to rod informations :" << endl;
  Info << "  Faces vertical to rod Number : " << face_vertical_to_rod_number_ << endl;
  Info << "  Total face area : " << total_face_area_vertical_to_rod_ << "  (m^2)" << nl << endl;


  forAll(faces_area_vector.boundaryField(), boundaryI)
  {
    if (faces_area_vector.boundaryField()[boundaryI].patch().name() == "INLET")
    {
      Info << " Boundary information: "
           << faces_area_vector.boundaryField()[boundaryI].patch().name() << endl;

      scalar total_face_area_ = 0.0;

      forAll(faces_area_mag.boundaryField()[boundaryI], faceI)
      {
        total_face_area_ += faces_area_mag.boundaryField()[boundaryI][faceI];
      }

      Info << "  Face Number       : " << faces_area_vector.boundaryField()[boundaryI].size() << endl;
      Info << "  Total Face Areas  : " << total_face_area_ << "  (m^2)" << nl << endl;
    }
  }


  forAll(faces_area_vector.boundaryField(), boundaryI)
  {
    if (faces_area_vector.boundaryField()[boundaryI].patch().name() == "OUTLET")
    {
      Info << " Boundary information: "
           << faces_area_vector.boundaryField()[boundaryI].patch().name() << endl;

      scalar total_face_area_ = 0.0;

      forAll(faces_area_mag.boundaryField()[boundaryI], faceI)
      {
        total_face_area_ += faces_area_mag.boundaryField()[boundaryI][faceI];
      }

      Info << "  Face Number       : " << faces_area_vector.boundaryField()[boundaryI].size() << endl;
      Info << "  Total Face Areas  : " << total_face_area_ << "  (m^2)" << nl << endl;
    }
  }


  if (rodsBuildOrNot)
  {
    labelList rodsSectionList_ = *rodEdgesListPtr_;

    Info << "  Total Rods sections :  " << rodsSectionList_.size() << endl;
  }
}


