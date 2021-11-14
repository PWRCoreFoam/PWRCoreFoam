/*----------------------------------------------------------------------------*\
| Class:                                                                       |
|     Foam::subChannelMesh                                                     |
\*----------------------------------------------------------------------------*/

#ifndef subChannelMesh_H
#define subChannelMesh_H

#include "fvCFD.H"
#include "fvMesh.H"
#include "IOdictionary.H"
#include "pointFieldFwd.H"


//****************************************************************************//

namespace Foam
{

/*----------------------------------------------------------------------------*\
|                       Class subChannelMesh Declaration                       |
\*----------------------------------------------------------------------------*/

class subChannelMesh
:
    public fvMesh
{

public:

    // Public typedefs


private:

    // Private data

        mutable labelList * rodEdgesListPtr_;
        mutable scalarList * rodEdgesLengthListPtr_;
        mutable pointField * rodEdgesCenterListPtr_;

        mutable labelList * rodEdgesListOnAssemblyEdgePtr_;
        mutable labelList * rodEdgesListOnCoreEdgePtr_;

        mutable labelList * gapListPtr_;

        label CoreMapDimension_;
        label AssemblyRodsDimension_;

        scalar AssemblyWidth_;

        label CenterAssemblyXIndex_;
        label CenterAssemblyYIndex_;

        scalar center_rod_X_index_;
        scalar center_rod_Y_index_;

        // Assembly index x and y
        labelListList * rod_assembly_index_X_Y_;

        // Rod index x and y in an assembly
        labelListList * rod_rod_index_X_Y_;

        // Core type for bundle or grid
        labelList * core_type_;

        labelList * cellRodsNumberList_;
        labelListList * cellRodsList_;

        labelListList core_assembly_distribution_;
        labelListList assembly_rod_layout_;

        scalarListList core_assembly_power_distribution_;
        scalarListList rod_pin_power_distribution_factor_;

        labelList core_type_z_direction_;
        scalarList core_height_z_direction_;

        scalar RodOuterDiameter_;

        scalar pitch_;
        Switch correct_rod_arrangement_;

        Switch test_case_with_assembly_box_;

        scalar assembly_box_corner_radius_;

        scalar mixing_coefficient_;

        bool subChannelMeshBuildOrNot;
        bool rodsBuildOrNot;
        bool mesh_modified_or_not;

        mutable scalarList flow_area_;
        mutable scalarList hydraulic_diameter_;


    // Private member functions

        void findRodEdgesList() const;
        void findRodEdgesOnAssemblyEdgeList() const;
        void findRodEdgesOnCoreEdgeList() const;
        void findGapList();
        void build_cellRodsList();
        void calc_index_assembly_and_rod();
        void calcRodEdgesLengthList();
        void calcRodEdgesCenterList();
        void calc_core_z_position();

        bool checkOnAssemblyEdge(vector p) const;

        bool check_face_vertical_to_z(vector p) const;

        scalar get_rod_section_length(label edgeI);
        vector get_rod_section_center(label edgeI);

        vector vector_remove_a_scalar(vector p, scalar s);

        void get_boundary_information_by_name(word  boundaryName,
                                              label & boundaryStart,
                                              label & boundarySize);

        label get_boundary_label_by_name(word boundaryName);

        void calc_flow_area();
        void calc_hydraulic_diameter();


public:

    // Runtime type information

        TypeName("subChannelMesh");


        // Construct from IOobject
        explicit subChannelMesh(const IOobject& io);

    // Destructor


    // Public member functions

        void test_function();

        void modify_mesh();
        void modify_mesh_for_core_edge();
        void modify_mesh_for_assembly_edge();

        void output_mesh_information();

        const labelList & get_rodEdgesList() const;
        const labelList & get_rodEdgesListOnAssemblyEdge() const;
        const labelList & get_rodEdgesListOnCoreEdge() const;
        const labelList & get_gapList() const;
        const scalarList & get_rodEdgesLengthList() const;
        const pointField & get_rodEdgesCenterList() const;

        inline const labelListList & get_rod_assembly_index_x_y() const;
        inline const labelListList & get_rod_rod_index_x_y() const;

        scalar get_rod_edge_z_max(label edge_i) const;
        scalar get_rod_edge_z_min(label edge_i) const;

        void readReactorCore();
        void buildSubChannelMesh();

        void createRods();

        inline label coreAssemblyDimension();
        inline label coreAssemblyDimension() const;

        inline label assemblyRodsDimension();
        inline label assemblyRodsDimension() const;

        inline scalar assemblyWidth();
        inline scalar assemblyWidth() const;

        inline label centerAssemblyXIndex();
        inline label centerAssemblyXIndex() const;

        inline label centerAssemblyYIndex();
        inline label centerAssemblyYIndex() const;

        inline label n_rods_sections();
        inline label n_rods_sections() const;

        inline scalarList & get_flow_area() const;
        inline scalarList & get_hydraulic_diameter() const;
        inline labelList & get_core_type() const;

        inline const scalarListList & get_core_assembly_power_distribution_factor() const;
        inline const scalarListList & get_rod_pin_power_distribution_factor() const;

        inline scalar mixing_coefficient();
        inline scalar mixing_coefficient() const;

        inline scalar get_pitch();
        inline scalar get_pitch() const;


    // Member operators


    // IOstream operators


};

}  // End namespace Foam

//****************************************************************************//

#include "subChannelMeshI.H"

//****************************************************************************//

#endif


