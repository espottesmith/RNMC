#pragma once
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <fstream>
#include <map>
#include <cassert> 
#include "stdio.h"
#include "stdlib.h"
#include <string>
#include "math.h"
#include <iostream>
#include <iostream>


#define DELTALOCAL 10000
#define DELTA 32768
#define EPSILON 0.0001

// universal defines inside namespace

#define FLERR __FILE__,__LINE__

#ifndef MIN
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#endif
#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

struct Site {

    Site() {}                               // default constructor does nothing

    Site(uint32_t i_in, uint32_t j_in,      // custom constructor for convenience 
        uint32_t k_in, float x_in, float y_in, float z_in,
        int species_in,
        bool can_adsorb_in) :
        i(i_in), j(j_in), k(k_in),
        x(x_in), y(y_in), z(z_in), 
        species(species_in), can_adsorb(can_adsorb_in) { }

    uint32_t i;                             // site location on lattice
    uint32_t j;
    uint32_t k;
    float x;                                // location in space
    float y;
    float z;
    int species;                            // species at site
    bool can_adsorb;                        // is the site in contact with the electrolyte?
    std::vector<int> shared_sites;           // indices of sites that also share this species

};

class Lattice {

private:

    /* ----------------------- handle memory ---------------------------------- */
    
    void *smalloc(int nbytes, const char *name);              // safe allocate         

    template <typename TYPE>
    TYPE *create(TYPE *&array, int n, const char *name);         // create 1D array

    template <typename TYPE>
    TYPE **create(TYPE **&array, int n1, int n2, const char *name); // create 2D array

    template <typename TYPE>
    void destroy(TYPE **array);                                  // destroy 2D array

    void sfree(void *ptr);                                        // safe free

    
    
    /* ----------------------- structural information ------------------------ */
    
    float latconst;                               // lattice constant
    float xlo, xhi, ylo,                   // bounding of box
          yhi, zlo, zhi;                       // (yscale * read in value)
    int ilo, ihi, klo, khi, jlo, jhi;                // geometry info neighbors
    bool is_xperiodic, is_yperiodic, is_zperiodic;          // 0 =   non-periodic, 1 = periodic
    
    int nsites;                                 // number of sites
    int nmax;                                   // max # sites, idneigh, numneigh can store at a given time

    int maxneigh;                               // max neighbors per site
    float maxz;                                   // the largest distance the lattice goes 

public: 

    std::unordered_map<int, Site> sites;                                // list of Sites for lattice
    std::unordered_map<int, uint32_t*> idneigh;                         // neighbor IDs for each site
    std::unordered_map<int, uint32_t> numneigh;                         // # of neighbors of each site
    std::unordered_map<int, char> edges;
    std::map<std::tuple<uint32_t, uint32_t, uint32_t>, int> loc_map;  // Mapping from site location (i,j,k) to site ID


    Lattice(float latconst_in, 
        int ilo_in, int ihi_in, int jlo_in,
        int jhi_in, int klo_in, int khi_in, 
        bool xperiodic_in, bool yperiodic_in, bool zperiodic_in);   
    
    Lattice(const Lattice& other);                          // copy constructor

    ~Lattice();

    void structured_lattice();
    
    void structured_connectivity();
    
    void offsets_3d(int **cmapone);
    
    void add_site(uint32_t i_in, uint32_t j_in, 
                  uint32_t k_in, float x_in, float y_in, float z_in,
                  bool can_adsorb_in, bool update_neighbors_in, bool meta_neighbors_in);

    void delete_site(int id);

    void update_neighbors(uint32_t n, bool meta_neighbors_in);

    float get_latconst();

    // TODO: make general for all types of periodicity 
    float get_maxz();

    // fill lattice with specified values from file
    void fill(std::string filename);

};

/* ---------------------------------------------------------------------- */

Lattice::Lattice(float latconst_in, 
        int ilo_in, int ihi_in, int jlo_in,
        int jhi_in, int klo_in, int khi_in, 
        bool is_xperiodic_in, bool is_yperiodic_in, bool is_zperiodic_in)  {

    latconst = latconst_in;
    
    // region of simulation input * lattice spacing
    xlo = ilo_in * latconst;
    xhi = ihi_in * latconst;
    ylo = jlo_in * latconst;
    yhi = jhi_in * latconst;
    zlo = klo_in * latconst;
    zhi = khi_in * latconst;
    
    // 0 = non-periodic, 1 = periodic
    is_xperiodic = is_xperiodic_in;
    is_yperiodic = is_yperiodic_in;
    is_zperiodic = is_zperiodic_in;
    
    nsites = 0;
    maxz = 0;
    nmax = DELTA;

    maxneigh = 6;

    // create sites on lattice
    structured_lattice();
    
    // set neighbors of each site
    structured_connectivity();
    
    /*for(int n = 0; n < nsites; ++n) {
        std::cout << "id: " << n << " [" << sites[n].x << ", " <<
        sites[n].y << ", " << sites[n].z << "]" << std::endl;
    }*/

} // Lattice()

/* ---------------------------------------------------------------------- */

Lattice::Lattice(const Lattice& other) {

    latconst = other.latconst;                               
    ilo = other.ilo;
    ihi = other.ihi;
    jlo = other.jlo;               
    jhi = other.jhi;
    klo = other.klo;
    khi = other.khi;                       
    xlo = other.xlo;
    xhi = other.xhi;
    ylo = other.ylo;
    yhi = other.yhi;
    zlo = other.zlo;
    zhi = other.zhi;                
    is_xperiodic = other.is_xperiodic;
    is_yperiodic = other.is_yperiodic;
    is_zperiodic = other.is_zperiodic;         
    
    for(auto it = other.sites.begin(); it != other.sites.end(); it++) {
        sites[it->first] = it->second;
    }

    nsites = other.nsites;                               
    nmax = other.nmax;                                

    maxneigh = other.maxneigh;                                             
        
    for(auto it = other.numneigh.begin(); it != other.numneigh.end(); it++) {
        numneigh[it->first] = it->second;
    }

    for(auto it = other.edges.begin(); it != other.edges.end(); it++) {
        edges[it->first] = it->second;
    }                       

    for(int i = 0; i < numneigh.size(); i++) {
        uint32_t* neighi;
        create(neighi, maxneigh, "create:neighi");
 
        for(size_t j = 0; j < other.numneigh.at(i); j++) {
            neighi[j] = other.idneigh.at(i)[j];
        }
        idneigh[i] = neighi;
    }                           


} // Lattice, custom constructor

/* ---------------------------------------------------------------------- */

Lattice::~Lattice() {

    for(size_t i = 0; i < idneigh.size(); i++) {
        sfree(idneigh[i]);
    }
} // ~Lattice()


/* ---------------------------------------------------------------------- */

void Lattice::structured_lattice() {
    
    // if not fully periodic IDs may be non-contiguous and/or ordered irregularly
    uint32_t nx, ny, nz;
    nx = (xhi - xlo / latconst);
    ny = (yhi - ylo / latconst);
    nz = (zhi - zlo / latconst);

    // if dim is periodic:
    //    lattice origin = lower box boundary
    //    loop bounds = 0 to N-1
    // if dim is non-periodic:
    //   lattice origin = 0.0
    //   loop bounds = enough to tile box completely, with all basis atoms
    
    if (is_xperiodic) {
        ilo = 0;
        ihi = nx-1;
      }
    else {
        ilo = (xlo / latconst);
        while ((ilo+1) * latconst > xlo) ilo--;
        ilo++;
        ihi = (xhi / latconst);
        while (ihi * latconst <= xhi) ihi++;
        ihi--;
      }

    if (is_yperiodic) {
        jlo = 0;
        jhi = ny-1;
    }
    else {
        jlo = (ylo / latconst);
        while ((jlo+1)*latconst > ylo) jlo--;
        jlo++;
        jhi = (yhi / latconst);
        while (jhi*latconst <= yhi) jhi++;
        jhi--;
    }

    if (is_zperiodic) {
        klo = 0;
        khi = nz-1;
    }
    else {
        klo = (zlo / latconst);
        while ((klo+1)*latconst > zlo) klo--;
        klo++;
        khi = (zhi / latconst);
        while (khi*latconst <= zhi) khi++;
        khi--;
    }

    
    // generate xyz coords and store them with site ID
    // tile the simulation box from origin, respecting PBC
    // site IDs should be contiguous if style = BOX and fully periodic
    // for non-periodic dims, check if site is within global box
    // for style = REGION, check if site is within region
    // if non-periodic or style = REGION, IDs may not be contiguous

    float x,y,z;
    bool can_adsorb;

    for (int k = klo; k <= khi; k++)
        for (int j = jlo; j <= jhi; j++)
            for (int i = ilo; i <= ihi; i++) {
                x = i*latconst;
                y = j*latconst;
                z = k*latconst;

                if (i == ihi && !is_xperiodic) {
                    can_adsorb = true;
                }
                else if (j == jhi && !is_yperiodic) {
                    can_adsorb = true;
                }
                else if (k == khi && !is_zperiodic) {
                    can_adsorb = true;
                }

                // By default, assume all lattice sites empty
                // TODO: This should use the global variable EMPTY_SITE
                // Don't update neighbors, since we'll use the connectivity function next
                add_site(i,j,k,x,y,z,can_adsorb,false,false);

    }


} // structered_lattice()

/* ----------------------------------------------------------------------
   generate site connectivity for on-lattice applications
 ------------------------------------------------------------------------- */

void Lattice::structured_connectivity() {

    int ineigh,jneigh,kneigh;
    uint32_t gid;
    int xneigh,yneigh,zneigh;
    
    float xprd = xhi - xlo;
    float yprd = yhi - ylo;
    float zprd = zhi - zlo;
    
    int nx = static_cast<int> (xprd / latconst);
    int ny = static_cast<int> (yprd / latconst);
    int nz = static_cast<int> (zprd / latconst);

    // create connectivity offsets
    
    int **cmap;                 // connectivity map for regular lattices
                                // cmap[maxneigh][3]
                                // 0,1,2 = i,j,k lattice unit cell offsets
    
    create(cmap,maxneigh,3,"create:cmap");

    offsets_3d(cmap);
    
    // generate global lattice connectivity for each site
    for (int i = 0; i < nsites; i++) {
        numneigh[i] = 0;

        uint32_t* neighi;
        create(neighi, maxneigh, "create:neighi");
        idneigh[i] = neighi;
      
        for (int neigh = 0; neigh < maxneigh; neigh++) {

            // ijkm neigh = indices of neighbor site
            // calculated from siteijk and cmap offsets

            ineigh = static_cast<int> (sites[i].x/latconst) + cmap[neigh][0];
            jneigh = static_cast<int> (sites[i].y/latconst) + cmap[neigh][1];
            kneigh = static_cast<int> (sites[i].z/latconst) + cmap[neigh][2];

            // xyz neigh = coords of neighbor site
            // calculated in same manner that structured_lattice() generated coords
            
            xneigh = static_cast<float> (ineigh) * latconst;
            yneigh = static_cast<float> (jneigh) * latconst;
            zneigh = static_cast<float> (kneigh) * latconst;
            
            // remap neighbor coords and indices into periodic box via ijk neigh
            // remap neighbor coords and indices into periodic box via ijk neigh

            if (is_xperiodic) {
                if (ineigh < 0) {
                    xneigh += xprd;
                    ineigh += nx;
                }
                if (ineigh >= nx) {
                    xneigh -= xprd;
                    xneigh = MAX(xneigh, static_cast<int> (xlo));
                    ineigh -= nx;
                }
            }
            if (is_yperiodic) {
                if (jneigh < 0) {
                    yneigh += yprd;
                    jneigh += ny;
                }
                if (jneigh >= ny) {
                    yneigh -= yprd;
                    yneigh = MAX(yneigh, static_cast<int> (ylo));
                    jneigh -= ny;
                }
            }
            if (is_zperiodic) {
                if (kneigh < 0) {
                    zneigh += zprd;
                    kneigh += nz;
                }
                if (kneigh >= nz) {
                    zneigh -= zprd;
                    zneigh = MAX(zneigh, static_cast<int> (zlo));
                    kneigh -= nz;
                }
            }

            // discard neighs that are outside non-periodic box or region
            if (!is_xperiodic && (xneigh < xlo || xneigh > xhi)) continue;
            if (!is_yperiodic && (yneigh < ylo || yneigh > yhi)) continue;
            if (!is_zperiodic && (zneigh < zlo || zneigh > zhi)) continue;

            // gid = global ID of neighbor
            // calculated in same manner that structured_lattice() generated IDs

            gid = ((kneigh-klo)*(jhi-jlo+1)*(ihi-ilo+1)) +
                  ((jneigh-jlo)*(ihi-ilo+1)) + ((ineigh-ilo));
            
            
        /*std::cout << "neighbor: (" << sites[gid].x << ", " << sites[gid].y << ", " 
                  << sites[gid].z << ") for: " << "[" << sites[i].x << ", " << 
                  sites[i].y << ", " << sites[i].z << "]" << std::endl;*/
            
        // add gid to neigh list of site i
        idneigh[i][numneigh[i]++] = gid;
      }
    }

   /* std::cout << "numneigh" << std::endl;
    for(int i = 0; i < nsites; i++) {
        std::cout << "[" << sites[i].x << ", " <<
        sites[i].y << ", " << sites[i].z << "]" << ",";
        std::cout << "num: " << numneigh[i] << std::endl;
    }
    
    std::cout << "idneigh" << std::endl;
    for(int i = 0; i < nsites; i++) {
        std::cout << "For site: " << i << "neighbors: ";
        for(int j = 0; j < numneigh[i]; j++) {
            std::cout << idneigh[i][j] << ", ";
        }
        std::cout << std::endl;
    }*/
    
    destroy(cmap);
} // structured_connectivity()

/* ---------------------------------------------------------------------- */

void Lattice::offsets_3d(int **cmapin) {
    
    int n = 0;
    double delx,dely,delz,r;
    double cutoff = latconst;

    for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
            for (int k = -1; k <= 1; k++) {
                delx = i * latconst;
                dely = j * latconst;
                delz = k * latconst;
                r = sqrt(delx * delx + dely * dely + delz * delz);
                if (r > cutoff - EPSILON && r < cutoff + EPSILON) {
                    assert(n != maxneigh);
                    cmapin[n][0] = i;
                    cmapin[n][1] = j;
                    cmapin[n][2] = k;
                    n++;
                }
            }
        }
    }
  assert(n == maxneigh);
} // offsets_3d

/* ---------------------------------------------------------------------- */

void Lattice::add_site(uint32_t i_in, uint32_t j_in, 
                       uint32_t k_in, float x_in, float y_in, float z_in,
                       bool can_adsorb_in, bool update_neighbors_in, bool meta_neighbors_in) {
    
    std::tuple<uint32_t, uint32_t, uint32_t> key = {i_in, j_in, k_in};
    if(loc_map.find(key) != loc_map.end()) {
        // site already exists
        return;
    }

    if (nsites == nmax) {
        nmax += DELTA;
        //can_adsorb.reserve(nmax);
        //can_desorb.reserve(nmax);
    }

    // Initialize neighbor information for this new site
    numneigh[nsites] = 0;

    uint32_t* neighi;
    create(neighi, maxneigh, "create:neighi");
    idneigh[nsites]= neighi;

    // initially empty site, species = 0
    sites[nsites] = Site{i_in, j_in, k_in, x_in, y_in, z_in, 0, can_adsorb_in};
    
    loc_map[key] = nsites;
    
    if(can_adsorb_in) {
        edges[nsites] = 'a';
    }

    if (update_neighbors_in) {
        update_neighbors(nsites, meta_neighbors_in);
    }

    if (sites[nsites].x < xlo) {
        xlo = sites[nsites].x;
        ilo = sites[nsites].i;
    }
    else if (sites[nsites].x > xhi) {
        xhi = sites[nsites].x;
        ihi = sites[nsites].i;
    }

    if (sites[nsites].y < ylo) {
        ylo = sites[nsites].y;
        jlo = sites[nsites].j;
    }
    else if (sites[nsites].y > yhi) {
        yhi = sites[nsites].y;
        jhi = sites[nsites].j;
    }

    if (sites[nsites].z < zlo) {
        zlo = sites[nsites].z;
        klo = sites[nsites].k;
    }
    else if (sites[nsites].z > zhi) {
        zhi = sites[nsites].z;
        khi = sites[nsites].k;
    }

    nsites++;

    // update running max distance // TODO: make general for all types of periodicity 
    if(z_in > maxz) {
        maxz = z_in;
    }

} // add_site()

/* ---------------------------------------------------------------------- */

void Lattice::delete_site(int id) {

    loc_map.erase({sites[id].i, sites[id].j, sites[id].k});

    assert(sites.find(id) != sites.end());

    // update neighbors
    update_neighbors(id, true);

    // delete from other hashes
    numneigh.erase(id);
    sfree(idneigh[id]);
    idneigh.erase(id);

    // remove from sites 
    sites.erase(id);
    

}

/* ---------------------------------------------------------------------- */

void Lattice::update_neighbors(uint32_t n, bool meta_neighbors_in) {

    float xprd = xhi - xlo;
    float yprd = yhi - ylo;
    float zprd = zhi - zlo;
    
    int nx = static_cast<int> (xprd / latconst);
    int ny = static_cast<int> (yprd / latconst);
    int nz = static_cast<int> (zprd / latconst);

    uint32_t left, right, backward, forward, down, up;
    
    left = sites[n].i - 1;
    right = sites[n].i + 1;
    if (is_xperiodic) {
        if (left < 0) {
            left += nx;
        }
        if (right >= nx) {
            right -= nx;
        }
    }
    
    backward = sites[n].j - 1;
    forward = sites[n].j + 1;
    if (is_yperiodic) {
        if (backward < 0) {
            backward += ny;
        }
        if (forward >= ny) {
            forward -= ny;
        }
    }

    down = sites[n].k - 1;
    up = sites[n].k + 1;
    if (is_zperiodic) {
        if (down < 0) {
            down += nz;
        }
        if (up >= nz) {
            up -= nz;
        }
    }

    std::vector<std::tuple<uint32_t,uint32_t,uint32_t>> ijk;
    ijk.resize(6);

    ijk[0] = {left, sites[n].j, sites[n].k};
    ijk[1] = {right, sites[n].j, sites[n].k};
    ijk[2] = {sites[n].i, backward, sites[n].k};
    ijk[3] = {sites[n].i, forward, sites[n].k};
    ijk[4] = {sites[n].i, sites[n].j, down};
    ijk[5] = {sites[n].i, sites[n].j, up};
    

    std::map<std::tuple<uint32_t,uint32_t,uint32_t>, int>::iterator it;

    // rest numneigh, and idneigh
    numneigh[n] = 0;
    uint32_t* neighi;
    create(neighi, maxneigh, "create:neighi");
    idneigh[n] = neighi;

    for (int q = 0; q < 6; q++) {
        auto it = loc_map.find(ijk[q]);
        if (it != loc_map.end()) {
            idneigh[n][numneigh[n]++] = it->second;
            
            if(meta_neighbors_in) {
                update_neighbors(it->second, false);
            }

        }
    }

    if(numneigh[n] == 0) {
        delete_site(n);
    }

} // update_neighbors()

/* ----------------------------------------------------------------------
   safe malloc 
------------------------------------------------------------------------- */

void *Lattice::smalloc(int nbytes, const char *name)
{
    if (nbytes == 0) return NULL;

    void *ptr = malloc(size_t(nbytes));
    if(ptr == NULL) {
        std::cout << name << '\n';
        assert(false);
    }               
    return ptr;
}
/* ----------------------------------------------------------------------
   create a 1d array 
------------------------------------------------------------------------- */

template <typename TYPE>
TYPE *Lattice::create(TYPE *&array, int n, const char *name)
{
    int nbytes = ((int) sizeof(TYPE)) * n;
    array = (TYPE *) smalloc(nbytes, name);
    return array;
}

/* ----------------------------------------------------------------------
   create a 2d array 
------------------------------------------------------------------------- */

template <typename TYPE>
TYPE **Lattice::create(TYPE **&array, int n1, int n2, const char *name) 
{
    int nbytes = ((int) sizeof(TYPE)) * n1*n2;
    TYPE *data = (TYPE *) smalloc(nbytes,name);
    nbytes = ((int) sizeof(TYPE *)) * n1;
    array = (TYPE **) smalloc(nbytes,name);
    
    int n = 0;
    for (int i = 0; i < n1; i++) {
    array[i] = &data[n];
    n += n2;
    }
    return array;
}

/* ----------------------------------------------------------------------
   destroy a 2d array 
------------------------------------------------------------------------- */

template <typename TYPE>
void Lattice::destroy(TYPE **array)
{
    if (array == NULL) return;
    sfree(array[0]);
    sfree(array);
}

/* ----------------------------------------------------------------------
   safe free 
------------------------------------------------------------------------- */

void Lattice::sfree(void *ptr)
{
  if (ptr == NULL) return;
  free(ptr);
}

/* ---------------------------------------------------------------------- */

float Lattice::get_latconst() {
    return latconst;
}

/* ---------------------------------------------------------------------- */

// TODO: make general for all types of periodicity 
float Lattice::get_maxz() {
    return maxz;
}

/* ---------------------------------------------------------------------- */
void Lattice::fill(std::string filename) {

    std::ifstream fin;
    fin.open(filename);

    if(!fin.is_open()) {
        std::cout << "Failed to open file: " << filename << "\n";
        exit(EXIT_FAILURE);
    }

    char type;
    char junk;
    double i_in, j_in, k_in;
    int species;

    fin >> type;

    if(type == 'L') {
        while(fin >> junk >> i_in >> junk >> j_in >> junk >> k_in >> junk >> species) {
            std::tuple<uint32_t, uint32_t, uint32_t> key = {i_in, j_in, k_in};
            sites[loc_map[key]].species = species;
        }
    }
    else if(type == 'A') {
        
        for(int k = klo; k <= khi; k++) {
            fin >> junk;
            for(int i = ilo; i <= ihi; i++) {
                for(int j = jlo; j <= jhi; j++) {
                    fin >> species;
                    std::tuple<uint32_t, uint32_t, uint32_t> key = {i, j, k}; 
                    sites[loc_map[key]].species = species;
                }
            }
            fin >> junk;
        }

    }
    else {
        std::cout << "Incorrect type of input: " << type << "\n";
        exit(EXIT_FAILURE);
    }

}


/* ---------------------------------------------------------------------- */


// TESTING //
/*
int main(int argc, char **argv) {
    
    Lattice *lattice = new Lattice(1, 0, 3, 0, 3, 0, 3, false, false, false);

    std::ofstream myfile;
    myfile.open ("test.txt");
    myfile << lattice->sites.size() << std::endl;

    for(auto it = lattice->sites.begin(); it != lattice->sites.end(); it++) {
        myfile << "[" << it->second.x << ", " <<
        it->second.y << ", " << it->second.z << "]" << ":";

        std::tuple<uint32_t, uint32_t, uint32_t> key = {it->second.x, it->second.y, it->second.z};
        int i = lattice->loc_map[key];  

        for(int j = 0; j < lattice->numneigh[i]; j++) {
            int id = lattice->idneigh[i][j];
            myfile << "(" << lattice->sites[id].x << ", " <<
            lattice->sites[id].y << ", " << lattice->sites[id].z << ")" << ",";
        }
        myfile << std::endl;
    }
    
    lattice->add_site(3, 3, 4, 3, 3, 4, true, true, true);
    lattice->add_site(2, 3, 4, 2, 3, 4, true, true, true);
    
    myfile << std::endl;
    myfile << "added ... " << std::endl;
    myfile << std::endl;
    myfile << lattice->sites.size() << std::endl;

    for(auto it = lattice->sites.begin(); it != lattice->sites.end(); it++) {
        myfile << "[" << it->second.x << ", " <<
        it->second.y << ", " << it->second.z << "]" << ":";

        std::tuple<uint32_t, uint32_t, uint32_t> key = {it->second.x, it->second.y, it->second.z};
        int i = lattice->loc_map[key];  

        for(int j = 0; j < lattice->numneigh[i]; j++) {
            int id = lattice->idneigh[i][j];
            myfile << "(" << lattice->sites[id].x << ", " <<
            lattice->sites[id].y << ", " << lattice->sites[id].z << ")" << ",";
        }
        myfile << std::endl;
    }
    
    std::tuple<uint32_t, uint32_t, uint32_t> key = {3, 3, 4};
    int id = lattice->loc_map[key];    
    lattice->delete_site(id);

    key = {2, 3, 4};
    id = lattice->loc_map[key];    
    lattice->delete_site(id);
    myfile << std::endl;
    myfile << "deleting ... " << std::endl;
    myfile << std::endl;

    myfile << lattice->sites.size() << std::endl;

    for(auto it = lattice->sites.begin(); it != lattice->sites.end(); it++) {
        myfile << "[" << it->second.x << ", " <<
        it->second.y << ", " << it->second.z << "]" << ":";

        std::tuple<uint32_t, uint32_t, uint32_t> key = {it->second.x, it->second.y, it->second.z};
        int i = lattice->loc_map[key];  

        for(int j = 0; j < lattice->numneigh[i]; j++) {
            int id = lattice->idneigh[i][j];
            myfile << "(" << lattice->sites[id].x << ", " <<
            lattice->sites[id].y << ", " << lattice->sites[id].z << ")" << ",";
        }
        myfile << std::endl;
    }

    myfile.close();
    delete lattice;

    
}*/