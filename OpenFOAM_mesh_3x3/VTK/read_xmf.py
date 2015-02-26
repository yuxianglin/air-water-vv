#! /usr/bin/env python
import numpy
def findGridElement(xmf,MeshTag='Spatial_Domain',id_in_collection=-1,verbose=0):
    """Try to find the element of the xml tree xmf that holds a uniform
    grid with the name given in MeshTag by searching through Temporal
    Grid Collections and Grid Collections.

    If MeshTag isn't found, uses the first entry in the Domain
    """
    Domain = xmf.getroot()[-1]
    GridCollection = None
    Grid = None
    for collection in Domain:
        if 'Name' in collection.attrib and MeshTag in collection.attrib['Name']:
            GridCollection = collection
            break
    if GridCollection == None:
        GridCollection = Domain[0]
    if verbose > 0:
        print "Trying GridCollection.tag= %s" % (GridCollection.tag)
    if GridCollection.attrib['GridType'] == 'Collection':
        Grid = GridCollection[-1]
    elif GridCollection.attrib['GridType'] == 'Uniform':
        Grid = GridCollection
    assert Grid.tag == 'Grid'
    assert Grid.attrib['GridType'] == 'Uniform'

    return Grid
def readMeshXdmf(xmf_archive_base,heavy_file_base,MeshTag="Spatial_Domain",hasHDF5=True,verbose=0):
    """
    start trying to read an xdmf archive with name xmf_archive_base.xmf
    assumes heavy_file_base.h5 has heavy data
    root Element is Xdmf
      last child of Xdmf which should be a Domain Element
         find child of Domain that is a Temporal Grid Collection with a name containing MeshTag, if None use first collection
            last child of Temporal Grid Collection should be a Uniform Grid at final time
               Attribute (usually 1) of child is  Topology  
                  set elementTopologyName to Type
                  if Type != Mixed
                    get text attribute and read this entry from  hdf5 file
                    set nNodes_element based on Type, nElements_global from leading dimension of elementNodesArray
                    create elementNodes_offset from Type and flatten elementNodesArray 
                  else  
                    get text attribute and read this entry from  hdf5 file to place in into xdmf_topology
                    generate elementNodesArray from xdmf_topology, calculating the number of elements using
                      walk through xdmf_topology   
               Attribute (usually 2) of child is Geometry  --> load data into nodeArray
                   set nNodes_global from nodeArray
               If has Attribute nodeMaterials read this from hdf file, else set to default of all zeros
               If has Attribute elementMaterialTypes, read this from hdf file, else set to default of all zeros                 

    returns a BasicMeshInfo object with the minimal information read
    """
    ###information about allowed Xdmf topologies
    #Xdmf cell type id to Name
    topologyid2name = {2:'Polyline',4:'Triangle',5:'Quadrilateral',6:'Tetrahedron',8:'Wedge',9:'Hexahedron',
                       112:'Mixed'} #Mixed isn't actually used 0x070
    #Topology name to number of local nodes
    topology2nodes = {'Polyline':2,'Triangle':3,'Quadrilateral':4,'Tetrahedron':4,'Wedge':6,'Hexahedron':8}

    #for output
    class BasicMeshInfo:
        def __init__(self):
            self.nNodes_global     = None
            self.nodeArray         = None
            self.nodeMaterialTypes = None
            self.nNodes_element    = None
            self.nElements_global  = None
            self.elementTopologyName = None
            self.elementNodesArray = None
            self.elementNodes_offset  = None
            self.elementMaterialTypes = None 
            self.nNodes_owned         = None
            self.nElements_owned      = None
    MeshInfo = BasicMeshInfo()
    try:
        from xml.etree import ElementTree as ET
        import tables
        xmf = ET.parse(xmf_archive_base+'.xmf')
        hdf5= tables.openFile(heavy_file_base+'.h5',mode="r") 
        assert hasHDF5
        Grid = findGridElement(xmf,MeshTag,id_in_collection=-1,verbose=verbose)

        #Geometry first
        Topology = None; Geometry  = None; NodeMaterials= None; ElementMaterials = None
        for i,leaf in enumerate(Grid):
            if verbose > 3:
                print "Grid leaf %d tag= %s " % (i,leaf.tag)
            if leaf.tag == 'Topology':
                Topology = Grid[i]
                if verbose > 3:
                    print "Topology found in leaf %d " % i
            elif leaf.tag == 'Geometry':
                Geometry = Grid[i]
                if verbose > 3:
                    print "Geometry found in leaf %d " % i
            elif leaf.tag == 'Attribute' and leaf.attrib['Name'] == 'nodeMaterialTypes':
                NodeMaterials = Grid[i]
                if verbose > 3:
                    print "NodeMaterials found in leaf %d " % i
            elif leaf.tag == 'Attribute' and leaf.attrib['Name'] == 'elementMaterialTypes':
                ElementMaterials = Grid[i]
                if verbose > 3:
                    print "ElementMaterials found in leaf %d " % i
        assert Geometry != None
        entry = Geometry[0].text.split(':')[-1]
        if verbose > 1:
            print "Reading nodeArray from %s " % entry
        MeshInfo.nodeArray = hdf5.getNode(entry).read()
        MeshInfo.nNodes_global = MeshInfo.nodeArray.shape[0]
        if NodeMaterials != None:
            entry = NodeMaterials[0].text.split(':')[-1]
            if verbose > 1:
                print "Reading nodeMaterialTypes from %s " % entry
            MeshInfo.nodeMaterialTypes = hdf5.getNode(entry).read()
        else:
            MeshInfo.nodeMaterialTypes = numpy.zeros((MeshInfo.nNodes_global,),'i')
        assert Topology != None
        
        if 'Type' in Topology.attrib:
            MeshInfo.elementTopologyName = Topology.attrib['Type']
        elif 'TopologyType' in Topology.attrib:
            MeshInfo.elementTopologyName = Topology.attrib['TopologyType']
        assert MeshInfo.elementTopologyName != None

        if verbose > 0:
            print "elementTopologyName= %s " % MeshInfo.elementTopologyName
        assert MeshInfo.elementTopologyName in topologyid2name.values()
        if MeshInfo.elementTopologyName != 'Mixed':
            MeshInfo.nNodes_element = topology2nodes[MeshInfo.elementTopologyName]
            entry = Topology[0].text.split(':')[-1]
            if verbose > 1:
                print "Reading  elementNodesArray from %s " % entry
            MeshInfo.elementNodesArray = hdf5.getNode(entry).read()
            assert MeshInfo.elementNodesArray.shape[1] == MeshInfo.nNodes_element
            MeshInfo.nElements_global = MeshInfo.elementNodesArray.shape[0]
            if verbose > 1:
                print "nElements_global,nNodes_element= (%d,%d) " % (MeshInfo.nElements_global,MeshInfo.nNodes_element)
            MeshInfo.elementNodes_offset = numpy.arange(MeshInfo.nElements_global*MeshInfo.nNodes_element+1,step=MeshInfo.nNodes_element,dtype='i')
        else:
            entry = Topology[0].text.split(':')[-1]
            if verbose > 1:
                print "Reading xdmf_topology from %s " % entry
            xdmf_topology = hdf5.getNode(entry).read()
            #build elementNodesArray and offsets now
            MeshInfo.nElements_global = 0
            i = 0
            while i < len(xdmf_topology):
                MeshInfo.nElements_global += 1
                nNodes_local = topology2nodes[topologyid2name[xdmf_topology[i]]]
                i += nNodes_local+1
            #
            if verbose > 2:
                print "Mixed topology found %s elements " % MeshInfo.nElements_global
            MeshInfo.elementNodes_offset = numpy.zeros((MeshInfo.nElements_global+1,),'i')

            i = 0; eN = 0
            while i < len(xdmf_topology):
                nNodes_local = topology2nodes[topologyid2name[xdmf_topology[i]]]
                MeshInfo.elementNodes_offset[eN+1] = MeshInfo.elementNodes_offset[eN] + nNodes_local
                eN += 1; i += nNodes_local+1
            MeshInfo.elementNodesArray = numpy.zeros((MeshInfo.elementNodes_offset[MeshInfo.nElements_global],),'i')
            i = 0; eN = 0
            while i < len(self.xdmf_topology):
                nNodes_local = topology2nodes[topologyid2name[xdmf_topology[i]]]
                MeshInfo.elementNodesArray[MeshInfo.elementNodes_offset[eN]:MeshInfo.elementNodes_offset[eN+1]][:] = xdmf_topology[i+1:i+1+nNodes_local][:]
                eN += 1; i += nNodes_local+1

        #
        if ElementMaterials != None:
            entry = ElementMaterials[0].text.split(':')[-1]
            if verbose > 1:
                print "Reading elementMaterialTypes from %s " % entry
            MeshInfo.elementMaterialTypes = hdf5.getNode(entry).read()
        else:
            MeshInfo.elementMaterialTypes = numpy.zeros((MeshInfo.nElements_global,),'i')
        #
        ###only serial for now
        MeshInfo.nNodes_owned = MeshInfo.nNodes_global
        MeshInfo.nElements_owned = MeshInfo.nElements_global
        hdf5.close()
    except:
        print "read from {0} failed ".format(xmf_archive_base)
        raise
    return MeshInfo
#

def writeHexMesh(mesh_info,hexfile_base,index_base=1):
    """
    Write a hex mesh in Ido's format with base numbering index_base
    HEX
    nNodes_global nElements_global
    x0 y0 z0
    x1 y1 z1
    ...
    xN yN zN
    [n0 n1 n2 n3 n4 n5 n6 n7 mat0]
    [n0 n1 n2 n3 n4 n5 n6 n7 mat1]
    """
    assert mesh_info.elementTopologyName=='Hexahedron'

    header="""HEX
{nNodes_global} {nElements_global}
""".format(nNodes_global=mesh_info.nNodes_global,nElements_global=mesh_info.nElements_global)

    with open(hexfile_base+'.mesh','w') as mout:
        mout.write(header)
        numpy.savetxt(mout,mesh_info.nodeArray)
        #format the elements, appending element material type
        elems_with_mat = numpy.append(mesh_info.elementNodesArray,mesh_info.elementMaterialTypes.reshape(mesh_info.nElements_global,1),axis=1)
        elems_with_mat[:,:-1] += index_base
        numpy.savetxt(mout,elems_with_mat,fmt='%d')

import numpy.testing as npt
from nose.tools import ok_ as ok
from nose.tools import eq_ as eq

def test_3x3_cube(verbose=0):
    """
    Read sample openfoam mesh from aggelos and check that the basic information is correct
    """
    mesh_info = readMeshXdmf('test_read_3x3','test_read_3x3',verbose=0)

    eq(mesh_info.nElements_global,27)
    eq(mesh_info.nNodes_global,64)
    
    eq(mesh_info.nElements_global,mesh_info.nElements_owned)
    eq(mesh_info.nNodes_global,mesh_info.nNodes_owned)

    eq(mesh_info.nodeArray.shape,(mesh_info.nNodes_owned,3))
    eq(mesh_info.elementNodesArray.shape,(mesh_info.nElements_owned,8))
    
    eq(mesh_info.elementTopologyName,'Hexahedron')

    eq(len(mesh_info.nodeMaterialTypes),mesh_info.nNodes_owned)
    eq(len(mesh_info.elementMaterialTypes),mesh_info.nElements_owned)

def test_write_3x3_cube(verbose=0):
    """
    Read sample openfoam mesh from aggelos and try to write in Ido's hex format
    """
    mesh_info = readMeshXdmf('test_read_3x3','test_read_3x3',verbose=0)

    writeHexMesh(mesh_info,'hexMesh_3x3',index_base=1)

if __name__ == '__main__':
    import nose
    nose.main(defaultTest='read_xmf:test_3x3_cube,read_xmf:test_write_3x3_cube')
