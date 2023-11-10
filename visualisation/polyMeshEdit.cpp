#include <iostream>
#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/helpers/Shortcuts.h>
#include <DGtal/helpers/ShortcutsGeometry.h>
#include "DGtal/io/readers/MeshReader.h"
#include "DGtal/io/writers/MeshWriter.h"

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/pick.h>
#include <polyscope/point_cloud.h>

#include "CLI11.hpp"


using namespace DGtal;
using namespace Z3i;

typedef PolygonalSurface<Z3i::RealPoint> PolySurface;

static PolySurface currentPolysurf;
static PolySurface firstPolysurf;

static std::vector<int> vectSelection;
static float paintRad = 1.0;
static float noiseLevel = 1.0;
static int partialF = 1;
static int randLarge = 100000;

static const int labelSelect = 200;
static const int labelRef = 50;

static std::string outputFileName {"result.obj"};


void addSurfaceInPolyscope(PolySurface &psurf){
    polyscope::removeStructure("InputMesh");
    std::vector<std::vector<std::size_t>> faces;
    vectSelection.clear();
    for(auto &face: psurf.allFaces()){
        faces.push_back(psurf.verticesAroundFace( face ));
        vectSelection.push_back(labelSelect);
    }
    
    auto digsurf = polyscope::registerSurfaceMesh("InputMesh",
                                                  psurf.positions(), faces);
    digsurf->addFaceScalarQuantity("selection", vectSelection);
    digsurf->setAllQuantitiesEnabled(true);
}

static Z3i::RealPoint
getFaceBarycenter(PolySurface &polysurff, const PolySurface::Face &aFace){
    Z3i::RealPoint res(0.0, 0.0, 0.0);
    for(auto const &v: polysurff.verticesAroundFace(aFace)){
        res += polysurff.positions()[v];
    }
    return res/polysurff.verticesAroundFace(aFace).size();
}

// Helper function
std::vector<PolySurface::Face> faceAround(PolySurface &polysurff,
                                          PolySurface::Face faceId,
                                          double radius){
    std::vector<PolySurface::Face> result;
    std::queue<PolySurface::Vertex> q;
    std::map<PolySurface::Face, bool> fVisited;
    std::map<PolySurface::Face, bool> vVisited;

    for (auto const &v : polysurff.verticesAroundFace(faceId) ){
        q.push(v);
    }
    fVisited[faceId] = true;
    bool addNewFaces = true;
    while (!q.empty() ) {
        PolySurface::Vertex v = q.front(); q.pop();
        auto listFace = polysurff.facesAroundVertex(v);
        for (auto const & f : listFace){
            if (fVisited.count(f)==0){
                if((getFaceBarycenter(polysurff, f) -
                    getFaceBarycenter(polysurff, faceId)).norm() < radius){
                    fVisited[f]=true;
                    result.push_back(f);
                    for (auto const & v : polysurff.verticesAroundFace(f)){
                        if (vVisited.count(v)==0){
                            vVisited[v]=true;
                            q.push(v);
                        }
                    }
                }
            }
        }
    }
    
    return result;
}
/**
 * Remove face from selectio with a probability of 1/selFreq
 */
void partialSelect(int selFreq=1){
    srand((unsigned) time(NULL));
    for ( unsigned int i = 0; i< currentPolysurf.nbFaces(); i++ ){
        if (vectSelection[i]==labelRef){
            if (rand()%selFreq==0){
                vectSelection[i]=labelSelect;
            }
        }
    }
}

void noisify(double scale = 0.01){
    srand((unsigned) time(NULL));
    for ( unsigned int i = 0; i< currentPolysurf.nbFaces(); i++ ){
        Z3i::RealPoint pDep((((double)(rand()%randLarge)-randLarge/2.0)/randLarge)*scale,
                            (((double)(rand()%randLarge)-randLarge/2.0)/randLarge)*scale,
                            (((double)(rand()%randLarge)-randLarge/2.0)/randLarge)*scale);
        if (vectSelection[i]==labelRef){
            for (auto f: currentPolysurf.verticesAroundFace(i)){
                currentPolysurf.positions()[f][0] += pDep[0];
                currentPolysurf.positions()[f][1] += pDep[1];
                currentPolysurf.positions()[f][2] += pDep[2];
            }
        }
    }
    addSurfaceInPolyscope(currentPolysurf);
}


DGtal::Mesh<DGtal::Z3i::RealPoint>
meshFixVertexUnity(DGtal::Mesh<DGtal::Z3i::RealPoint> &aMesh){
    DGtal::Mesh<DGtal::Z3i::RealPoint> res;
    std::map<DGtal::Z3i::RealPoint, bool> mapPt;
    std::vector<bool> vertexUsed (aMesh.nbVertex(), false);
    for ( unsigned int f = 0; f< aMesh.nbFaces(); f++ ){
        auto face = aMesh.getFace(f);
        for (unsigned int i = 0; i<face.size(); i++){
            vertexUsed[face[i]] = true;
            if (mapPt.count(aMesh.getVertex(face[i]))==0){
                mapPt[aMesh.getVertex(face[i])]=true;
            }
        }
    }
    
    std::vector<int> translateIndexId;
    int currentIndex = 0;
    for(unsigned int i = 0; i < aMesh.nbVertex(); i++ ){
        if (vertexUsed[i]){
            res.addVertex(aMesh.getVertex(i));
            translateIndexId.push_back(currentIndex);
            currentIndex++;
        }else{
            translateIndexId.push_back(-1);
        }
    }
    vertexUsed =  std::vector<bool> (aMesh.nbVertex(), false);
    
    for ( unsigned int f = 0; f< aMesh.nbFaces(); f++ ){
        auto face = aMesh.getFace(f);
        for (unsigned int i = 0; i<face.size(); i++){
            face[i]=translateIndexId[face[i]];
            vertexUsed[face[i]] = true;
        }
        res.addFace(face);
    }
    return res;
}
            

            
void deleteSelectedFaces(){
    PolySurface newSur;
    std::vector<bool> vertexUsed (currentPolysurf.nbVertices(), false);
    for ( unsigned int i = 0; i< currentPolysurf.nbFaces(); i++ ){
        if (vectSelection[i]==labelSelect){
            auto face = currentPolysurf.verticesAroundFace( i );
            for (auto v: face){
                vertexUsed[v] = true;
            }
        }
    }
    auto lp = currentPolysurf.positions();
    for(unsigned int i = 0; i < currentPolysurf.nbVertices(); i++){
        if(vertexUsed[i]){
            newSur.addVertex(lp[i]);
        }
    }
    std::vector<int> translateIndexId;
    int currentIndex = 0;
    for(unsigned int i = 0; i < currentPolysurf.nbVertices(); i++ ){
        if (vertexUsed[i]){
            translateIndexId.push_back(currentIndex);
            currentIndex++;
        }else{
            translateIndexId.push_back(-1);
        }
    }
    for ( unsigned int f = 0; f< currentPolysurf.nbFaces(); f++ ){
        if (vectSelection[f]==labelSelect){
            auto face = currentPolysurf.verticesAroundFace( f );
            for (unsigned int i = 0; i<face.size(); i++){
                face[i]=translateIndexId[face[i]];
            }
            newSur.addPolygonalFace(PolySurface::PolygonalFace( face.cbegin(), face.cend()));
        }
    }
    newSur.build();
    currentPolysurf = newSur;
    addSurfaceInPolyscope(newSur);
}

void callbackFaceID() {
    srand((unsigned) time(NULL));

    ImGui::Text("Setting selection size:");
    ImGui::SliderFloat("radius values", &paintRad, 0.01f, 10.0f, "size = %.3f");

    ImGui::Separator();

    ImGui::Text("Set selection freq:");
    ImGui::SliderInt(" freq (1=select all, 2=select 1over2)", &partialF, 1, 10, "freq = %i");
    ImGui::Separator();
    ImGui::Text("Noise parameters:");
    ImGui::SliderFloat("noise scale", &noiseLevel, 0.001, 2, "scale = %f");
    ImGui::Separator();

    ImGui::Text("Action:");
    if (ImGui::Button("Clear selection")) {
        for(auto &i : vectSelection){
            i=labelSelect;
        }
    }
    ImGui::SameLine();
    if (ImGui::Button("delete selected faces")) {
        deleteSelectedFaces();
    }
    ImGui::SameLine();
    if (ImGui::Button("noisify selected faces")) {
        noisify(noiseLevel);
    }
    ImGui::Separator();
    ImGui::Text("IO");

    if (ImGui::Button("save in .obj")) {
        std::ofstream obj_stream( outputFileName.c_str() );
        MeshHelpers::exportOBJ(obj_stream, currentPolysurf);
    }
    ImGui::SameLine();

    if (ImGui::Button("reload src")) {
        currentPolysurf = firstPolysurf;
        addSurfaceInPolyscope(currentPolysurf);
        
    }
    ImGui::Separator();
    ImGui::Text("Polyscope interface:");

    if (ImGui::Button("show ")) {
        polyscope::options::buildGui=true;
    }
    ImGui::SameLine();
    if (ImGui::Button("hide")) {
        polyscope::options::buildGui=false;
    }

    ImGuiIO& io = ImGui::GetIO();
    auto digsurf = polyscope::getSurfaceMesh("InputMesh");
    if (io.MouseDoubleClicked[0]) {
        unsigned long indexSelect = polyscope::pick::getSelection().second;
        unsigned long nb = 0;
        // face selected
        if (indexSelect >= currentPolysurf.nbVertices()){
            nb = (unsigned long) polyscope::pick::getSelection().second - currentPolysurf.nbVertices();
        }else{
            // vertex selected (selecting a face connected to it)
            if(currentPolysurf.facesAroundVertex(polyscope::pick::getSelection().second)
               .size()> 0){
                nb = currentPolysurf.facesAroundVertex(polyscope::pick::getSelection().second)[0];
            }
        }
        
        if (nb > 0 && nb < vectSelection.size()){
            auto fVois = faceAround(currentPolysurf, nb, paintRad);
            vectSelection[nb] = 0;
            srand((unsigned) time(NULL));
            for (auto f: fVois){
                if (rand()%partialF==0){
                    vectSelection[f]=labelRef;
                }else{
                    vectSelection[f]=labelSelect;
                }
             }
        }
    }
    digsurf -> removeQuantity("selection");
    digsurf -> addFaceScalarQuantity("selection", vectSelection);
}





int main(int argc, char** argv)
{
    std::string inputFileName {""};
    std::string colorFileName; // The file containing the index to be colored in the colors
    
    // parse command line using CLI ----------------------------------------------
    CLI::App app;
    app.description("polyMeshEdit \n"
                    " polyMeshEdit -i file.obj  \n");
    app.add_option("-i,--input,1", inputFileName, "an input mesh file in .obj or .off format." )
    ->required()
    ->check(CLI::ExistingFile);
    app.add_option("-o,--output,2", outputFileName, "an output mesh file in .obj or .off format." )
    ->check(CLI::ExistingFile);
    
    
    app.get_formatter()->column_width(40);
    CLI11_PARSE(app, argc, argv);
    polyscope::init();
    polyscope::options::programName = "PolyMeshEdit";
    polyscope::options::buildGui=false;
    // read input mesh
    DGtal::Mesh<DGtal::Z3i::RealPoint> aMesh(true);
    aMesh << inputFileName;
    aMesh.removeIsolatedVertices();
    DGtal::MeshHelpers::mesh2PolygonalSurface(aMesh, currentPolysurf);
    polyscope::state::userCallback = callbackFaceID;
    addSurfaceInPolyscope(currentPolysurf);
    firstPolysurf = currentPolysurf;
    polyscope::show();
    return 0;
    
}







