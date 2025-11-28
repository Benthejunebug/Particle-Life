#include "Universe.h"
#include "HSV.h"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <thread>
#include <future>
#include <functional>
#include <math.h>
#include <stdlib.h>
#include <mutex>

static const float RADIUS = 6.0f;
static const float DIAMETER = 2.0f * RADIUS;
static const float R_SMOOTH = 2.0f;
static const int   NEARBYRESETCOUNT = 1000;

struct arg_struct {
    size_t start;
    size_t end;
};

Universe::Universe(size_t num_types, size_t num_particles, int width, int height) {
    //Initialize everything
    
    //Old seeds
    //1089903674
    m_rand_gen.seed(57454860);
    //m_rand_gen.seed((unsigned int)time(0));
    SetPopulation(num_types, num_particles);
    SetSize(float(width), float(height));
    //numthreads = std::thread::hardware_concurrency();
    numthreads = 1;
    m_center_x = m_width * 0.5f;
    m_center_y = m_height * 0.5f;
    m_zoom = 1.0f;
    m_attract_mean = 0.0f;
    m_attract_std = 0.0f;
    m_minr_lower = 0.0f;
    m_minr_upper = 0.0f;
    m_maxr_lower = 0.0f;
    m_maxr_upper = 0.0f;
    m_friction = 0.0f;
    m_flat_force = false;
    m_wrap = false;
    nearbyCounter = 0;
    dualExplosions = false;
    singleExplosions = false;
    randomDecay = false;
    randomDecayFactor = 100;
    th.reserve(numthreads);
    srand((unsigned int)time(0));
    SetRandomTypeAllocation();
    
    for (auto i = randomTypeAssignment.begin(); i !=  randomTypeAssignment.end(); ++i){
        //std::cout << *i << std::endl;
    }
    bonding = false;
    bondChart = {};
    bondTable = {};
    
    //Player particle implementation
    plyrP.x = m_center_x;
    plyrP.y = m_center_y;
    plyrP.vx = 0;
    plyrP.vy = 0;
    plyrP.e_orbitals = 1;
    plyrDirection = 1;
    plyrBonding = false;
    plyrBondsMax = 20;
}

// Goal: to rewrite this as memory safe just getting rid of the locks
// do I need to make struts?


void Universe::addBond(size_t p1I, size_t p2I){
    
    for(auto bond: bondChart){
        if((bond[0] == p1I && bond[1] == p2I) ||
           (bond[1] == p1I && bond[0] == p2I)){
            return;
        }
    }
    bondChart.push_back({p1I,p2I});
    m_particles[p1I].e_orbitals += m_types.ENumber(m_particles[p2I].type);
    m_particles[p2I].e_orbitals += m_types.ENumber(m_particles[p1I].type);
}


void Universe::removeBond(size_t p1I, size_t p2I){
    int i = 0;
    m_particles[p1I].e_orbitals -= m_types.ENumber(m_particles[p2I].type);
    m_particles[p2I].e_orbitals -= m_types.ENumber(m_particles[p1I].type);
    for (auto bond: bondChart){
        if((bond[0] == p1I && bond[1] == p2I) ||
           (bond[1] == p1I && bond[0] == p2I)){
            bondChart.erase(bondChart.begin()+i);
            return;
        }
        i++;
    }
    return;
}

void Universe::addBondTwo(size_t p1I, size_t p2I){
    std::lock_guard<std::mutex> guard(bond_mutex_p[p1I]);
    std::lock_guard<std::mutex> guard2(bond_mutex_p[p2I]);
    
    m_particles[p1I].bonds.push_back(p2I);
    m_particles[p2I].bonds.push_back(p1I);
    m_particles[p1I].e_orbitals += m_types.ENumber(m_particles[p2I].type);
    m_particles[p2I].e_orbitals += m_types.ENumber(m_particles[p1I].type);
}

void Universe::removeBondTwo(size_t p1I, size_t p2I){
    std::lock_guard<std::mutex> guard(bond_mutex_p[p1I]);
    std::lock_guard<std::mutex> guard2(bond_mutex_p[p2I]);
    
    Particle& p = m_particles[p1I];
    Particle& q = m_particles[p2I];
    p.e_orbitals -= m_types.ENumber(q.type);
    q.e_orbitals -= m_types.ENumber(p.type);
    
    int i = 0;
    for (auto bond: p.bonds){
        if(bond == p2I){
            p.bonds.erase(p.bonds.begin()+i);
            break;
        }
        i++;
    }
    
    i = 0;
    for (auto bond: q.bonds){
        if(bond == p1I){
            q.bonds.erase(q.bonds.begin()+i);
            return;
        }
        i++;
    }
}

//void Universe::addBondThree(size_t p1I, size_t p2I){
//    std::lock_guard<std::mutex> guard(bond_mutex_p[p1I]);
//    std::lock_guard<std::mutex> guard2(bond_mutex_p[p2I]);
//    
//    bondTable[p1I].row[p2I] = true;
//    bondTable[p2I].row[p1I] = true;
//    m_particles[p1I].e_orbitals += m_types.ENumber(m_particles[p2I].type);
//    m_particles[p2I].e_orbitals += m_types.ENumber(m_particles[p1I].type);
//}
//
//void Universe::removeBondThree(size_t p1I, size_t p2I){
//    std::lock_guard<std::mutex> guard(bond_mutex_p[p1I]);
//    std::lock_guard<std::mutex> guard2(bond_mutex_p[p2I]);
//    
//    Particle& p = m_particles[p1I];
//    Particle& q = m_particles[p2I];
//    p.e_orbitals -= m_types.ENumber(q.type);
//    q.e_orbitals -= m_types.ENumber(p.type);
//    
//    bondTable[p1I].row[p2I] = false;
//    bondTable[p2I].row[p1I] = false;
//}

void Universe::addPlyrBond(size_t p){
    plyrBondList.push_back(p);
    m_particles[p].e_orbitals += 1;
}

void Universe::removePlyrBond(size_t p){
    int i = 0;
    m_particles[p].e_orbitals -= 1;
    for ( size_t bond: plyrBondList){
        if(bond == p){
            plyrBondList.erase(plyrBondList.begin()+i);
            return;
        }
        i++;
    }
}


bool Universe::bonded(size_t particleIndex){
    std::lock_guard<std::mutex> guard(bond_mutex);
    
    for (auto bond: bondChart){
        if (bond[0] == particleIndex){
            return true;
        } else if (bond[1] == particleIndex){
            return true;
        }
    }
    return false;
}


bool Universe::isPlyrBonded(size_t p){
    for (int i=0;i<plyrBondList.size();i++){
        if(plyrBondList[i] == p){
            return true;
        }
    }
    return false;
}

bool isPlyrBonded(std::vector<size_t> plyrBondList,size_t p){
    for (int i=0;i<plyrBondList.size();i++){
        if(plyrBondList[i] == p){
            return true;
        }
    }
    return false;
}

//std::vector<bool> Universe::getBondPartnersThree(size_t particleIndex){
//    bond_mutex_p[particleIndex].lock();
//    BondRow bondPartners = bondTable[particleIndex];
//    bond_mutex_p[particleIndex].unlock();
//    return bondPartners.row;
//}

std::vector<size_t> Universe::getBondPartners(size_t particleIndex){
    std::lock_guard<std::mutex> guard(bond_mutex);
    std::vector<size_t> bondPartners = {};
    for (int i=0;i<bondChart.size();i++){
        if (bondChart[i][0] == particleIndex){
            bondPartners.push_back(bondChart[i][1]);
        }else if(bondChart[i][1] == particleIndex){
            bondPartners.push_back(bondChart[i][0]);
        }
    }
    return bondPartners;
}

std::vector<size_t> getBondPartners(std::vector<std::vector<size_t>> bondChart, size_t particleIndex){
    std::vector<size_t> bondPartners = {};
    for (int i=0;i<bondChart.size();i++){
        if (bondChart[i][0] == particleIndex){
            bondPartners.push_back(bondChart[i][1]);
        }else if(bondChart[i][1] == particleIndex){
            bondPartners.push_back(bondChart[i][0]);
        }
    }
    return bondPartners;
}


int Universe::getBondPartnersCount(size_t particleIndex){
    std::lock_guard<std::mutex> guard(bond_mutex);
    int bondPartners = 0;
    for (int i=0;i<bondChart.size();i++){
        if (bondChart[i][0] == particleIndex){
            bondPartners += 1;
        }else if(bondChart[i][1] == particleIndex){
            bondPartners += 1;
        }
    }
    return bondPartners;
}

bool Universe::CheckTypeAssignment(std::vector<int> randomTypes){
    //std::cout << "#############" << std::endl;
    for (int i=0; i<randomTypes.size(); i++){
        //std::cout << randomTypes.at(i) << std::endl;
        
        if (i == randomTypes.at(i) || i == randomTypes.at(randomTypes.at(i))){
            //std::cout << randomTypes.at(i) << std::endl;
            //std::cout << "FALSE" << std::endl;
            return false;
        }
    }
    return true;
}

std::vector<int> Universe::RandomizeTypes(int size){
    std::vector<int> shuffledTypes = {};
    for (int i=0;i<size;i++){
        shuffledTypes.push_back(i);
    }
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(shuffledTypes.begin(), shuffledTypes.end(), g);
    for (auto i = shuffledTypes.begin(); i !=  shuffledTypes.end(); ++i){
        //std::cout << *i << std::endl;
    }
    
    return shuffledTypes;
}

void Universe::ReSeed(float attract_mean, float attract_std, float minr_lower, float minr_upper,
                      float maxr_lower, float maxr_upper, float friction, bool flat_force) {
    m_attract_mean = attract_mean;
    m_attract_std = attract_std;
    m_minr_lower = minr_lower;
    m_minr_upper = minr_upper;
    m_maxr_lower = maxr_lower;
    m_maxr_upper = maxr_upper;
    m_friction = friction;
    m_flat_force = flat_force;
    SetRandomTypes();
    SetRandomParticles();
    SetRandomTypeAllocation();
    std::lock_guard<std::mutex> guard(bond_mutex);
    bondChart.clear();
    bondTable.clear();
    plyrBondList.clear();
    
    
    //Player Particle Implementation
    plyrRadius = 0.8f*minr_upper;
    plyrMinR = plyrRadius+5;
    plyrMaxR = maxr_upper;
    plyrAttractionMagnitude = abs(attract_mean);
    plyrAttractionIncrement = attract_std/5;
    
}

void Universe::SetRandomTypeAllocation(){
    randomTypeAssignment = RandomizeTypes((int)m_types.Size());
    while (!CheckTypeAssignment(randomTypeAssignment)){
        randomTypeAssignment = RandomizeTypes((int)m_types.Size());
    }
}

void Universe::SetPopulation(size_t num_types, size_t num_particles) {
    m_types.Resize(num_types);
    m_particles.resize(num_particles);
    bond_mutex_p = new std::mutex[(int)num_particles];
    bondTable = std::vector<BondRow>(num_particles, BondRow(num_particles));
}

void Universe::SetRandomTypes() {
    std::normal_distribution<float>       rand_attr(m_attract_mean, m_attract_std);
    std::uniform_real_distribution<float> rand_minr(m_minr_lower, m_minr_upper);
    std::uniform_real_distribution<float> rand_maxr(m_maxr_lower, m_maxr_upper);
    
    std::cout << "-----------------------" << std::endl;
    
    for (size_t i = 0; i < m_types.Size(); ++i) {
        m_types.Color(i) = FromHSV(float(i) / m_types.Size(), 1.0f, float(i % 2)*0.25f + 0.75f);
        m_types.ENumber(i) = rand()%8+1;
        std::cout << "ENUMMM: "<< m_types.ENumber(i) << std::endl;
        for (size_t j = 0; j < m_types.Size(); ++j) {
            if (i == j) {
                m_types.Attaract(i, j) = -std::abs(rand_attr(m_rand_gen));
                m_types.MinR(i, j) = DIAMETER;
            } else {
                m_types.Attaract(i, j) = rand_attr(m_rand_gen);
                m_types.MinR(i, j) = std::max(rand_minr(m_rand_gen), DIAMETER);
            }
            m_types.MaxR(i, j) = std::max(rand_maxr(m_rand_gen), m_types.MinR(i, j));
            
            //Keep radii symmetric
            m_types.MaxR(j, i) = m_types.MaxR(i, j);
            m_types.MinR(j, i) = m_types.MinR(i, j);
        }
    }
}

void Universe::SetRandomParticles() {
    std::uniform_int_distribution<int>    rand_type(0, int(m_types.Size() - 1));
    std::uniform_real_distribution<float> rand_uni(0.0f, 1.0f);
    std::normal_distribution<float>       rand_norm(0.0f, 1.0f);
    
    //Player implementation
    plyrV = 0.5f;
    
    for (size_t i = 0; i < m_particles.size(); ++i) {
        Particle& p = m_particles[i];
        p.type = uint8_t(rand_type(m_rand_gen));
        p.x = (rand_uni(m_rand_gen)*0.06f + 0.5f) * m_width;
        p.y = (rand_uni(m_rand_gen)*0.06f + 0.5f) * m_height;
        p.vx = rand_norm(m_rand_gen) * 0.2f;
        p.vy = rand_norm(m_rand_gen) * 0.2f;
        p.e_orbitals = m_types.ENumber(p.type);
    }
}


void Universe::chunkInteractionsTwo(size_t start, size_t end) {
    for (size_t i = start; i < end; ++i) {
        //Current particle
        Particle& p = m_particles[i];
        
        bond_mutex_p[i].lock();
        std::vector<size_t> bonds = p.bonds;
        bond_mutex_p[i].unlock();
        
        
        //Interactions
        //for (size_t j = 0; j < p.nearbyParticles.size(); ++j) {
        unsigned long int last = m_particles.size();
        if(playerParticleActivated){
            last += 1;
        }
        for (size_t j = 0; j < last; ++j) {
            if(i==j){
                continue;
            }
            bool isPlyrP = (j == m_particles.size())? true : false;
            
            //Other particle
            Particle& q = isPlyrP? plyrP : m_particles[j];
            
            
            bool bonded = false;
            if (!isPlyrP){
                for (auto bond: bonds){
                    if (bond == j){
                        bonded = true;
                    }
                }
            } else {
                for(auto bond: plyrBondList){
                    if (bond ==  i){
                        bonded = true;
                    }
                }
            }
            
            //Get deltas
            float dx = q.x - p.x;
            float dy = q.y - p.y;
            if (m_wrap) {
                if (dx > m_center_x) {
                    dx -= m_width;
                } else if (dx < -m_center_x) {
                    dx += m_width;
                }
                if (dy > m_center_y) {
                    dy -= m_height;
                } else if (dy < -m_center_y) {
                    dy += m_height;
                }
            }
            
            float minR = m_types.MinR(p.type, q.type);
            float maxR = m_types.MaxR(p.type, q.type);
            
            float maxR2 = maxR*maxR;
            
            if ((std::max(dx, dy) > maxR2) && !bonded) {
                continue;
            }
            
            //Get distance squared
            const float r2 = dx*dx + dy*dy;
            
            
         
            
            if (isPlyrP){
                minR = plyrMinR;
                maxR = plyrMaxR;
            }
            
            if ((r2 > maxR*maxR || r2 < 0.01f) && !bonded) {
                continue;
            }
            
            float bondFactor = (abs(m_types.Attaract(p.type, q.type)+m_attract_mean)/(m_attract_mean+1.5*m_attract_std))/2;
            //            if (bondFactor > 0.5f){
            //                bondFactor /= bondFactor*1.5f;
            //            }
            
            int qENum = m_types.ENumber(q.type);
            
            if(isPlyrP){
                bondFactor = (60+(rand()%40))/100.0f;
                qENum = 1;
            }
            
            int peorbsAdjust = bonded? p.e_orbitals-qENum : p.e_orbitals;
            float bondFactor2 = ((abs(peorbsAdjust+qENum)%9)/8.0f);
            
            const float bondR = 1.02*bondFactor*minR*bondFactor2; //minR/100;
            
            //            if ( rand()%100000000 < 2){
            //                std::cout << "++++++++++++++"<< std::endl;
            //                std::cout << "Electron Orbital: " << p.e_orbitals << std::endl;
            //                std::cout << "p factor: "<< abs(p.e_orbitals-8)/8.0f << std::endl;
            //                std::cout << "q factor: " <<abs(q.e_orbitals-8)/8.0f << std::endl;
            //                std::cout << "BONDS: " << bonds.size() << std::endl;
            //                std::cout << "BOND R: " << bondR << std::endl;
            //                std::cout << "MIN R: " << minR << std::endl;
            //               std::cout << "MAX R: " << maxR << std::endl;
            //                std::cout << "attract: " << m_types.Attaract(p.type, q.type)<< std::endl;
            //
            //               std::cout << "bond factor 2: " <<bondFactor2<< std::endl;
            //
            //
            //           }
            //
          
            
            //Normalize displacement
            const float r = std::sqrt(r2);
            dx /= r;
            dy /= r;
            
            bool brokeBonds = false;
            
            if (bonding && !isPlyrP){
                //Make Bond
                int p1eN = m_types.ENumber(p.type);
                int p2eN = m_types.ENumber(q.type);
                if(p.e_orbitals+p2eN <= 8 && q.e_orbitals+p1eN <= 8){
                    if (!bonded &&  r < bondR){
                        bonded = true;
                        addBondTwo(i, j);
                    }
                }
                
            }
            
            if (bonding && isPlyrP && plyrBonding && plyrBondList.size() < plyrBondsMax){
                if (!bonded &&  r < bondR){
                    bonded = true;
                    addPlyrBond(i);
                }
            }
            
            //Break Bond
            if(!bonded && r<bondR/10 && rand()<1000){
                brokeBonds = true;
                bond_mutex_p[i].lock();
                for(size_t bond: bonds){
                    bond_mutex_p[bond].lock();
                    std::vector<size_t>& otherBonds = m_particles[bond].bonds;
                    int l = 0;
                    for(size_t ob: otherBonds){
                        if(ob == i){
                            otherBonds.erase(otherBonds.begin()+l);
                            break;
                        }
                        l++;
                    }
                    bond_mutex_p[bond].unlock();
                }
                if (isPlyrBonded(i)){
                    removePlyrBond(i);
                }
                p.bonds.empty();
                bond_mutex_p[i].unlock();
            }
            
            if(bonded){
                
                if((r>maxR*(bondFactor2) && (isPlyrP? r>30+plyrMaxR/((((int)i*3)/m_particles.size())+1) : true)) || p.e_orbitals>8 || q.e_orbitals>8){
                    brokeBonds = true;
                    if(isPlyrP){
                        removePlyrBond(i);
                    }else{
                        removeBondTwo(i, j);
                    }
                    bonded = false;
                }
            }
            
            
            //Calculate force
            
            float f = 0.0f;
            if (r > minR) {
                if (m_flat_force) {
                    f = m_types.Attaract(p.type, q.type);
                } else {
                    const float numer = 2.0f * std::abs(r - 0.5f*(maxR + minR));
                    const float denom = maxR - minR;
                    float attract = isPlyrP? plyrAttractDeNone*plyrAttractionMagnitude : m_types.Attaract(p.type, q.type);
                    f =  attract * (1.0f - numer / denom);
                }
                
            } else {
                f = R_SMOOTH*minR*(1.0f/(minR + R_SMOOTH) - 1.0f / (r + R_SMOOTH));
            }
            
            float bondPoint = 15.0f;
            
            if(isPlyrP){
                bondPoint = plyrMaxR/((((int)i*3)/m_particles.size())+1.5);
            }
            
            if (bonding && bonded){
                //if (r-bondPoint < 0){
                //f = 0.9f*(1/(1/(bondPoint+2)-1/(r+2)));//
                //f = 100*pow(bondPoint-r, 3);
                f = 0.09f*bondFactor2*(r-bondPoint);
            }
            
            
            /////////////////////////////////////////////////////////////
            if (!isPlyrP){
                if (dualExplosions && !bonded){
                    if (r < minR/2.5 && rand()%100>95 /*&& p.type == q.type*/){
                        //std::cout << "start type:" << (int)p.type << std::endl;
                        std::uniform_int_distribution<int>    rand_type(0, int(m_types.Size() - 1));
                        if(p.type != q.type){
                            p.type = uint8_t((int)(p.type + ((p.type+q.type)+(maxR/minR)*100))%m_types.Size());
                        }
                        //p.type = (uint8_t)rand_type(m_rand_gen);
                        
                        //std::cout << "end type:" << (int)p.type << std::endl;
                        
                    }
                }
                
                if (singleExplosions && !bonded){
                    if (r<minR/2 && rand()%100 > 98){
                        if ((int)p.type < m_types.Size()-1){
                            p.type += 1;
                        }else{
                            p.type = uint8_t(0);
                        }
                    } else if (r<minR/2.6 && rand()%1000 > 9995){
                        p.type = (uint8_t)((int)(p.type+pow(maxR/minR, 2))%m_types.Size());
                    }
                }
                
                if (randomDecay){
                    if (rand()<randomDecayFactor){
                        if ((int)p.type < m_types.Size()-1){
                            try{
                                p.type = randomTypeAssignment[p.type];
                            }catch(std::exception e){
                                std::cout << "FAIL: " << p.type << std::endl;
                                std::cout << (int)p.type << std::endl;
                            }
                        }else{
                            p.type = uint8_t(0);
                        }
                    }
                }
            }
            if(brokeBonds){
                f*=3;
            }
            //Apply force
            p.vx += f * dx;
            p.vy += f * dy;
        }
        
        
        
        
    }
    
    
}
//
//void Universe::chunkInteractionsThree(size_t start, size_t end) {
//    for (size_t i = start; i < end; ++i) {
//        //Current particle
//        Particle& p = m_particles[i];
//        std::vector<bool> bonds = getBondPartnersThree(i);
//        
//        
//        //Interactions
//        //for (size_t j = 0; j < p.nearbyParticles.size(); ++j) {
//        unsigned long int last = m_particles.size();
//        if(playerParticleActivated){
//            last += 1;
//        }
//        
//        for (size_t j = 0; j < last; ++j) {
//            if(i==j){
//                continue;
//            }
//            bool isPlyrP = (j == m_particles.size())? true : false;
//            
//            //Other particle
//            Particle& q = isPlyrP? plyrP : m_particles[j];
//            
//            
//            bool bonded = false;
//            if (!isPlyrP){
//                bonded = bonds[j];
//            } else {
//                for(auto bond: plyrBondList){
//                    if (bond ==  i){
//                        bonded = true;
//                    }
//                }
//            }
//            
//            //Get deltas
//            float dx = q.x - p.x;
//            float dy = q.y - p.y;
//            if (m_wrap) {
//                if (dx > m_center_x) {
//                    dx -= m_width;
//                } else if (dx < -m_center_x) {
//                    dx += m_width;
//                }
//                if (dy > m_center_y) {
//                    dy -= m_height;
//                } else if (dy < -m_center_y) {
//                    dy += m_height;
//                }
//            }
//            
//            //Get distance squared
//            const float r2 = dx*dx + dy*dy;
//            
//            
//            float minR = m_types.MinR(p.type, q.type);
//            float maxR = m_types.MaxR(p.type, q.type);
//            
//            if (isPlyrP){
//                minR = plyrMinR;
//                maxR = plyrMaxR;
//            }
//            
//            if ((r2 > maxR*maxR || r2 < 0.01f) && !bonded) {
//                continue;
//            }
//            
//            float bondFactor = (abs(m_types.Attaract(p.type, q.type)+m_attract_mean)/(m_attract_mean+1.5*m_attract_std))/2;
//            //            if (bondFactor > 0.5f){
//            //                bondFactor /= bondFactor*1.5f;
//            //            }
//            
//            int qENum = m_types.ENumber(q.type);
//            
//            if(isPlyrP){
//                bondFactor = (60+(rand()%40))/100.0f;
//                qENum = 1;
//            }
//            
//            int peorbsAdjust = bonded? p.e_orbitals-qENum : p.e_orbitals;
//            float bondFactor2 = ((abs(peorbsAdjust+qENum)%9)/8.0f);
//            
//            const float bondR = 1.02*bondFactor*minR*bondFactor2; //minR/100;
//            
//            //            if ( rand()%100000000 < 2){
//            //                std::cout << "++++++++++++++"<< std::endl;
//            //                std::cout << "Electron Orbital: " << p.e_orbitals << std::endl;
//            //                std::cout << "p factor: "<< abs(p.e_orbitals-8)/8.0f << std::endl;
//            //                std::cout << "q factor: " <<abs(q.e_orbitals-8)/8.0f << std::endl;
//            //                std::cout << "BONDS: " << bonds.size() << std::endl;
//            //                std::cout << "BOND R: " << bondR << std::endl;
//            //                std::cout << "MIN R: " << minR << std::endl;
//            //               std::cout << "MAX R: " << maxR << std::endl;
//            //                std::cout << "attract: " << m_types.Attaract(p.type, q.type)<< std::endl;
//            //
//            //               std::cout << "bond factor 2: " <<bondFactor2<< std::endl;
//            //
//            //
//            //           }
//            //
//            
//            
//            //Normalize displacement
//            const float r = std::sqrt(r2);
//            dx /= r;
//            dy /= r;
//            
//            bool brokeBonds = false;
//            
//            if (bonding && !isPlyrP){
//                //Make Bond
//                int p1eN = m_types.ENumber(p.type);
//                int p2eN = m_types.ENumber(q.type);
//                if(p.e_orbitals+p2eN <= 8 && q.e_orbitals+p1eN <= 8){
//                    if (!bonded &&  r < bondR){
//                        bonded = true;
//                        addBondThree(i, j);
//                    }
//                }
//                
//            }
//            
//            if (bonding && isPlyrP && plyrBonding && plyrBondList.size() < plyrBondsMax){
//                if (!bonded &&  r < bondR){
//                    bonded = true;
//                    addPlyrBond(i);
//                }
//            }
//            
//            //Break Bond
//            if(!bonded && r<bondR/10 && rand()<1000){
//                brokeBonds = true;
//                bond_mutex_p[i].lock();
//                for(int k = 0;k<bonds.size();k++){
//                    bond_mutex_p[k].lock();
//                    bondTable[k].row[i];
//                    bond_mutex_p[k].unlock();
//                }
//                if (isPlyrBonded(i)){
//                    removePlyrBond(i);
//                }
//                bondTable[i] = BondRow(m_particles.size());
//                bond_mutex_p[i].unlock();
//            }
//            
//            if(bonded){
//                
//                if((r>maxR*(bondFactor2) && (isPlyrP? r>30+plyrMaxR/((((int)i*3)/m_particles.size())+1) : true)) || p.e_orbitals>8 || q.e_orbitals>8){
//                    brokeBonds = true;
//                    if(isPlyrP){
//                        removePlyrBond(i);
//                    }else{
//                        removeBondThree(i, j);
//                    }
//                    bonded = false;
//                }
//            }
//            
//            
//            //Calculate force
//            
//            float f = 0.0f;
//            if (r > minR) {
//                if (m_flat_force) {
//                    f = m_types.Attaract(p.type, q.type);
//                } else {
//                    const float numer = 2.0f * std::abs(r - 0.5f*(maxR + minR));
//                    const float denom = maxR - minR;
//                    float attract = isPlyrP? plyrAttractDeNone*plyrAttractionMagnitude : m_types.Attaract(p.type, q.type);
//                    f =  attract * (1.0f - numer / denom);
//                }
//                
//            } else {
//                f = R_SMOOTH*minR*(1.0f/(minR + R_SMOOTH) - 1.0f / (r + R_SMOOTH));
//            }
//            
//            float bondPoint = 15.0f;
//            
//            if(isPlyrP){
//                bondPoint = plyrMaxR/((((int)i*3)/m_particles.size())+1.5);
//            }
//            
//            if (bonding && bonded){
//                //if (r-bondPoint < 0){
//                //f = 0.9f*(1/(1/(bondPoint+2)-1/(r+2)));//
//                //f = 100*pow(bondPoint-r, 3);
//                f = 0.09f*bondFactor2*(r-bondPoint);
//            }
//            
//            
//            /////////////////////////////////////////////////////////////
//            if (!isPlyrP){
//                if (dualExplosions && !bonded){
//                    if (r < minR/2.5 && rand()%100>95 /*&& p.type == q.type*/){
//                        //std::cout << "start type:" << (int)p.type << std::endl;
//                        std::uniform_int_distribution<int>    rand_type(0, int(m_types.Size() - 1));
//                        if(p.type != q.type){
//                            p.type = uint8_t((int)(p.type + ((p.type+q.type)+(maxR/minR)*100))%m_types.Size());
//                        }
//                        //p.type = (uint8_t)rand_type(m_rand_gen);
//                        
//                        //std::cout << "end type:" << (int)p.type << std::endl;
//                        
//                    }
//                }
//                
//                if (singleExplosions && !bonded){
//                    if (r<minR/2 && rand()%100 > 98){
//                        if ((int)p.type < m_types.Size()-1){
//                            p.type += 1;
//                        }else{
//                            p.type = uint8_t(0);
//                        }
//                    } else if (r<minR/2.6 && rand()%1000 > 9995){
//                        p.type = (uint8_t)((int)(p.type+pow(maxR/minR, 2))%m_types.Size());
//                    }
//                }
//                
//                if (randomDecay){
//                    if (rand()<randomDecayFactor){
//                        if ((int)p.type < m_types.Size()-1){
//                            try{
//                                p.type = randomTypeAssignment[p.type];
//                            }catch(std::exception e){
//                                std::cout << "FAIL: " << p.type << std::endl;
//                                std::cout << (int)p.type << std::endl;
//                            }
//                        }else{
//                            p.type = uint8_t(0);
//                        }
//                    }
//                }
//            }
//            if(brokeBonds){
//                f*=3;
//            }
//            //Apply force
//            p.vx += f * dx;
//            p.vy += f * dy;
//        }
//        
//        
//        
//        
//    }
//    
//    
//}
//

void Universe::PrintTypeCount(){
    std::cout << "###############" << std::endl;
    for (int i = 0; i<m_types.Size();i++){
        int count = 0;
        for (int j = 0; j<m_particles.size();j++){
            if ((int)m_particles[j].type == i){
                count++;
            }
        }
        std::cout << "Type " << i << ": " << count << std::endl;
    }
}


void Universe::updateParticle(Particle& p){
    
    //Update position and velocity
    p.x += p.vx;
    p.y += p.vy;
    p.vx *= (1.0f - m_friction);
    p.vy *= (1.0f - m_friction);
    
    //Check for wall collisions
    if (m_wrap) {
        if (p.x < 0) {
            p.x += m_width;
        } else if (p.x >= m_width) {
            p.x -= m_width;
        }
        if (p.y < 0) {
            p.y += m_height;
        } else if (p.y >= m_height) {
            p.y -= m_height;
        }
    } else {
        if (p.x <= DIAMETER) {
            p.vx = -p.vx;
            p.x = DIAMETER;
        } else if (p.x >= m_width - DIAMETER) {
            p.vx = -p.vx;
            p.x = m_width - DIAMETER;
        }
        if (p.y <= DIAMETER) {
            p.vy = -p.vy;
            p.y = DIAMETER;
        } else if (p.y >= m_height - DIAMETER) {
            p.vy = -p.vy;
            p.y = m_height - DIAMETER;
        }
    }
}

void Universe::chunkUpdate(size_t start, size_t end) {
    //Update position
    for (size_t i = start; i < end; ++i) {
        //Current particle
        Particle& p = m_particles[i];
        updateParticle(p);
    }
    
}

void Universe::Step() {
    //if (nearbyCounter == 0){ AssignNearbyParticles();}
    numthreads = std::thread::hardware_concurrency();
    th.reserve(numthreads);
    if (nearbyCounter<NEARBYRESETCOUNT){
        nearbyCounter+=1;
    }else{
        //PrintTypeCount();
        std::cout << "0000000000000000000000000000" << std::endl;
        std::cout << "BOND CHART: " << bondChart.size() << " bonds"<< std::endl;
        //        for (auto bond: bondChart){
        //            std::cout << bond[0] << ", " << bond[1] << std::endl;
        //        }
        
        nearbyCounter = 0;
    }
    
    //th.clear();
    
    size_t start = 0;
    
    size_t end = m_particles.size();
    size_t chunk = (numthreads>0 ? ((end-start)/numthreads) : end - start);
//    std::cout << chunk << std::endl;
    //std::cout << numthreads << std::endl;
    //chunkInteractions(start, end);
    std::vector<std::future<std::vector<std::tuple<size_t, size_t, bool>>>> resChunk = {};
    
    if(bondingMethod == 0){
        for( int i = 0; i<numthreads; i++){
            size_t s = start + i*chunk;
            size_t e = s + chunk;
            
            auto bondCopy(bondChart);
            
            resChunk.push_back(std::async(std::launch::async,&Universe::chunkInteractions,this, s, e, bondCopy, plyrBondList)); //{&Universe::chunkInteractions, this, s, e});
        }
        
        std::vector<std::tuple<size_t, size_t, bool>> changes = {};
        for(auto &f: resChunk){
            auto newchanges = f.get();
            changes.insert(changes.end(), newchanges.begin(), newchanges.end() );
//            for (auto change : newchanges){
//                std::cout << std::get<0>(change) << ", "
//                          << std::get<1>(change) << ", "
//                          << std::get<2>(change) << std::endl;
//            }
            
            //So it works when printing all that stuff. I thought at first that it maybe was b/c it took some time, but actually it may lazily compute this putting off joining the thread until the actual values are needed.
        }
        resChunk.clear();
        
        for(auto change : changes){
            bool isPlyrP = (std::get<0>(change) == m_particles.size())? true : false;
            if (std::get<2>(change)){
                isPlyrP? addPlyrBond(std::get<1>(change)) : addBond(std::get<0>(change), std::get<1>(change));
            }else {
                isPlyrP? removePlyrBond(std::get<1>(change)) : removeBond(std::get<0>(change), std::get<1>(change));
            }
            
        }
    }else if(bondingMethod == 1){
        for( int i = 0; i<numthreads; i++){
            size_t s = start + i*chunk;
            size_t e = s + chunk;
            th.push_back(std::thread{&Universe::chunkInteractionsTwo, this, s, e});
        }
        
        
        for(auto &t : th){
            t.join();
        }
        
        th.clear();

    }
//    } else {
//        for( int i = 0; i<numthreads; i++){
//            size_t s = start + i*chunk;
//            size_t e = s + chunk;
//            th.push_back(std::thread{&Universe::chunkInteractionsThree, this, s, e});
//        }
//    }
//    
 
    
    
    

    
//    // Update particle bond boolean
//    for(size_t i=0; i<m_particles.size(); i++){
//        Particle& p = m_particles[i];
//        if(newBondingCode){
//            if(p.bonds.size() > 0){
//                p.bonded = true;
//            }else{
//                p.bonded = false;
//            }
//        }else{
//            p.bonded = bonded(i);
//        }
//    }
    
    if (playerParticleActivated){
        //std::cout << SameDirection(plyrDirection, PLYRUP) << std::endl;
        float v = plyrV;
        //if two directions
        if (plyrDirection>7){
            v = sqrt(v);
        }
        
        if (SameDirection(plyrDirection, PLYRUP)){
            plyrP.vy += -v;
        } else  if (SameDirection(plyrDirection, PLYRDN)){
            plyrP.vy += v;
        }
        
        if (SameDirection(plyrDirection, PLYRLF)){
            plyrP.vx += -v;
        } else  if (SameDirection(plyrDirection, PLYRRT)){
            plyrP.vx += v;
        }
        
    }
    //ClearPlayerDirection();
    //update position
    for( int i = 0; i<numthreads;i++){
        size_t s = start + i*chunk;
        size_t e = s + chunk;
        th.push_back(std::thread{&Universe::chunkUpdate, this, s, e});
    }
    
    
    if(playerParticleActivated){
        updateParticle(plyrP);
    }
    
    for(auto &t : th){
        t.join();
    }
    th.clear();
    
    
    
    //        for (size_t i = 0; i < m_particles.size(); ++i) {
    //            //Current particle
    //            Particle& p = m_particles[i];
    //
    //            //Interactions
    //            for (size_t j = 0; j < m_particles.size(); ++j) {
    //                //Other particle
    //                const Particle& q = m_particles[j];
    //
    //                //Get deltas
    //                float dx = q.x - p.x;
    //                float dy = q.y - p.y;
    //                if (m_wrap) {
    //                    if (dx > m_width*0.5f) {
    //                        dx -= m_width;
    //                    } else if (dx < -m_width*0.5f) {
    //                        dx += m_width;
    //                    }
    //                    if (dy > m_height*0.5f) {
    //                        dy -= m_height;
    //                    } else if (dy < -m_height*0.5f) {
    //                        dy += m_height;
    //                    }
    //                }
    //
    //                //Get distance squared
    //                const float r2 = dx*dx + dy*dy;
    //                const float minR = m_types.MinR(p.type, q.type);
    //                const float maxR = m_types.MaxR(p.type, q.type);
    //
    //                if (r2 > maxR*maxR || r2 < 0.01f) {
    //                    continue;
    //                }
    //
    //                //Normalize displacement
    //                const float r = std::sqrt(r2);
    //                dx /= r;
    //                dy /= r;
    //
    //                //Calculate force
    //                float f = 0.0f;
    //                if (r > minR) {
    //                    if (m_flat_force) {
    //                        f = m_types.Attaract(p.type, q.type);
    //                    } else {
    //                        const float numer = 2.0f * std::abs(r - 0.5f*(maxR + minR));
    //                        const float denom = maxR - minR;
    //                        f = m_types.Attaract(p.type, q.type) * (1.0f - numer / denom);
    //                    }
    //                } else {
    //                    f = R_SMOOTH*minR*(1.0f/(minR + R_SMOOTH) - 1.0f / (r + R_SMOOTH));
    //                }
    //
    //                //Apply force
    //                p.vx += f * dx;
    //                p.vy += f * dy;
    //            }
    //        }
    //  //Update position
    //  for (size_t i = 0; i < m_particles.size(); ++i) {
    //    //Current particle
    //    Particle& p = m_particles[i];
    //
    //    //Update position and velocity
    //    p.x += p.vx;
    //    p.y += p.vy;
    //    p.vx *= (1.0f - m_friction);
    //    p.vy *= (1.0f - m_friction);
    //
    //    //Check for wall collisions
    //    if (m_wrap) {
    //      if (p.x < 0) {
    //        p.x += m_width;
    //      } else if (p.x >= m_width) {
    //        p.x -= m_width;
    //      }
    //      if (p.y < 0) {
    //        p.y += m_height;
    //      } else if (p.y >= m_height) {
    //        p.y -= m_height;
    //      }
    //    } else {
    //      if (p.x <= DIAMETER) {
    //        p.vx = -p.vx;
    //        p.x = DIAMETER;
    //      } else if (p.x >= m_width - DIAMETER) {
    //        p.vx = -p.vx;
    //        p.x = m_width - DIAMETER;
    //      }
    //      if (p.y <= DIAMETER) {
    //        p.vy = -p.vy;
    //        p.y = DIAMETER;
    //      } else if (p.y >= m_height - DIAMETER) {
    //        p.vy = -p.vy;
    //        p.y = m_height - DIAMETER;
    //      }
    //    }
    //  }
}

void Universe::Draw(sf::RenderWindow& window, float opacity) const {
    if(bondingMethod == 0){
        for(auto bond: bondChart){
            Particle p1 = m_particles[bond[0]];
            Particle p2 = m_particles[bond[1]];
            const float p1x = (p1.x - m_center_x)*m_zoom + float(m_width/2);
            const float p1y = (p1.y - m_center_y)*m_zoom + float(m_height/2);
            const float p2x = (p2.x - m_center_x)*m_zoom + float(m_width/2);
            const float p2y = (p2.y - m_center_y)*m_zoom + float(m_height/2);
            if ((p1x > m_center_x+m_width/4 && p2x < m_center_x-+m_width/4)||(p1x < m_center_x-m_width/4 && p2x > m_center_x+m_width/4) ||
                (p1y > m_center_y+m_height/4 && p2y < m_center_y-m_height/4)||(p1y < m_center_y-m_height/4 && p2y > m_center_y+m_height/4)){
                continue;
            }
            sf::VertexArray lines(sf::LinesStrip, 2);
            lines[0].position = sf::Vector2f(p1x, p1y);
            lines[0].color = sf::Color::White;
            lines[1].position = sf::Vector2f(p2x, p2y);
            lines[1].color = sf::Color::White;
            
            window.draw(lines);
        }
    }else if(bondingMethod == 1){
        for(Particle p: m_particles){
            for(size_t bond: p.bonds){
                Particle p1 = p;
                Particle p2 = m_particles[bond];
                const float p1x = (p1.x - m_center_x)*m_zoom + float(m_width/2);
                const float p1y = (p1.y - m_center_y)*m_zoom + float(m_height/2);
                const float p2x = (p2.x - m_center_x)*m_zoom + float(m_width/2);
                const float p2y = (p2.y - m_center_y)*m_zoom + float(m_height/2);
                if ((p1x > m_center_x+m_width/4 && p2x < m_center_x-+m_width/4)||(p1x < m_center_x-m_width/4 && p2x > m_center_x+m_width/4) ||
                    (p1y > m_center_y+m_height/4 && p2y < m_center_y-m_height/4)||(p1y < m_center_y-m_height/4 && p2y > m_center_y+m_height/4)){
                    continue;
                }
                sf::VertexArray lines(sf::LinesStrip, 2);
                lines[0].position = sf::Vector2f(p1x, p1y);
                lines[0].color = sf::Color::White;
                lines[1].position = sf::Vector2f(p2x, p2y);
                lines[1].color = sf::Color::White;
                
                window.draw(lines);
            }
        }
    } else {
        for(int i=0;i<m_particles.size(); i++){
            for(int j=0;j<m_particles.size(); j++){
                if (bondTable[i].row[j]){
                    Particle p1 = m_particles[i];
                    Particle p2 = m_particles[j];
                    const float p1x = (p1.x - m_center_x)*m_zoom + float(m_width/2);
                    const float p1y = (p1.y - m_center_y)*m_zoom + float(m_height/2);
                    const float p2x = (p2.x - m_center_x)*m_zoom + float(m_width/2);
                    const float p2y = (p2.y - m_center_y)*m_zoom + float(m_height/2);
                    if ((p1x > m_center_x+m_width/4 && p2x < m_center_x-+m_width/4)||(p1x < m_center_x-m_width/4 && p2x > m_center_x+m_width/4) ||
                        (p1y > m_center_y+m_height/4 && p2y < m_center_y-m_height/4)||(p1y < m_center_y-m_height/4 && p2y > m_center_y+m_height/4)){
                        continue;
                    }
                    sf::VertexArray lines(sf::LinesStrip, 2);
                    lines[0].position = sf::Vector2f(p1x, p1y);
                    lines[0].color = sf::Color::White;
                    lines[1].position = sf::Vector2f(p2x, p2y);
                    lines[1].color = sf::Color::White;
                    
                    window.draw(lines);
                }
            }
        }
    }
    
    
    for(auto bond: plyrBondList){
        Particle p1 = plyrP;
        Particle p2 = m_particles[bond];
        const float p1x = (p1.x - m_center_x)*m_zoom + float(m_width/2);
        const float p1y = (p1.y - m_center_y)*m_zoom + float(m_height/2);
        const float p2x = (p2.x - m_center_x)*m_zoom + float(m_width/2);
        const float p2y = (p2.y - m_center_y)*m_zoom + float(m_height/2);
        if ((p1x > m_center_x+m_width/4 && p2x < m_center_x-+m_width/4)||(p1x < m_center_x-m_width/4 && p2x > m_center_x+m_width/4) ||
            (p1y > m_center_y+m_height/4 && p2y < m_center_y-m_height/4)||(p1y < m_center_y-m_height/4 && p2y > m_center_y+m_height/4)){
            continue;
        }
        sf::VertexArray lines(sf::LinesStrip, 2);
        lines[0].position = sf::Vector2f(p1x, p1y);
        lines[0].color = sf::Color::White;
        lines[1].position = sf::Vector2f(p2x, p2y);
        lines[1].color = sf::Color::White;
        
        window.draw(lines);
    }
    
    sf::CircleShape circle;
    circle.setRadius(RADIUS * m_zoom);
    circle.setOrigin(circle.getRadius(), circle.getRadius());
    for (size_t i = 0; i < m_particles.size(); ++i) {
        const Particle& p = m_particles[i];
        const float x = (p.x - m_center_x)*m_zoom + float(m_width/2);
        const float y = (p.y - m_center_y)*m_zoom + float(m_height/2);
        circle.setPosition(x, y);
        sf::Color col = m_types.Color(p.type);
        col.a = uint8_t(opacity * 255);
        circle.setFillColor(col);
        window.draw(circle);
    }
    
    //Draw Player Particle
    if (playerParticleActivated){
        
        circle.setRadius(plyrRadius*m_zoom);
        circle.setOrigin(circle.getRadius(), circle.getRadius());
        const float x = (plyrP.x - m_center_x)*m_zoom + float(m_width/2);
        const float y = (plyrP.y - m_center_y)*m_zoom + float(m_height/2);
        circle.setPosition(x, y);
        circle.setFillColor(plyrColor);
        window.draw(circle);
    }
    
}


int Universe::GetIndex(int x, int y) const {
    float cx, cy;
    ToCenter(x, y, cx, cy);
    for (size_t i = 0; i < m_particles.size(); ++i) {
        const float dx = m_particles[i].x - cx;
        const float dy = m_particles[i].y - cy;
        if (dx*dx + dy*dy < RADIUS*RADIUS) {
            return int(i);
        }
    }
    return -1;
}

float Universe::GetParticleX(int index) const {
    return m_particles[index].x;
}

float Universe::GetParticleY(int index) const {
    return m_particles[index].y;
}

void Universe::ToCenter(int x, int y, float& cx, float& cy) const {
    cx = m_center_x + float(x - m_width / 2) / m_zoom;
    cy = m_center_y + float(y - m_height / 2) / m_zoom;
}

void Universe::Zoom(float cx, float cy, float zoom) {
    //Apply the zoom
    m_center_x = cx;
    m_center_y = cy;
    m_zoom = std::max(1.0f, zoom);
    
    //Clamp to make sure camera doesn't go out of bounds
    m_center_x = std::min(m_center_x, float(m_width) * (1.0f - 0.5f/m_zoom));
    m_center_y = std::min(m_center_y, float(m_height) * (1.0f - 0.5f / m_zoom));
    m_center_x = std::max(m_center_x, float(m_width) * (0.5f / m_zoom));
    m_center_y = std::max(m_center_y, float(m_height) * (0.5f / m_zoom));
}

void Universe::PrintParams() const {
    std::cout << "\nAttract:\n";
    for (size_t i = 0; i < m_types.Size(); ++i) {
        for (size_t j = 0; j < m_types.Size(); ++j) {
            std::cout << std::fixed << std::setw(8) << std::setprecision(4) << m_types.Attaract(i, j) << "  ";
        }
        std::cout << "\n";
    }
    std::cout << "MinR:\n";
    for (size_t i = 0; i < m_types.Size(); ++i) {
        for (size_t j = 0; j < m_types.Size(); ++j) {
            std::cout << std::fixed << std::setw(8) << std::setprecision(4) << m_types.MinR(i, j) << "  ";
        }
        std::cout << "\n";
    }
    std::cout << "MaxR:\n";
    for (size_t i = 0; i < m_types.Size(); ++i) {
        for (size_t j = 0; j < m_types.Size(); ++j) {
            std::cout << std::fixed << std::setw(8) << std::setprecision(4) << m_types.MaxR(i, j) << "  ";
        }
        std::cout << "\n";
    }
}

//CHUNK UPDATE BACKUP


std::vector<std::tuple<size_t, size_t, bool>>  Universe::chunkInteractions(size_t start, size_t end, std::vector<std::vector<size_t>> bondChart, std::vector<size_t> plyrBondList) {
    // true for add false for remove and if the value is at mparticles size
    // then it is our player particle
    std::vector<std::tuple<size_t, size_t, bool>> bondChanges = {};
    //std::cout << std::this_thread::get_id() << std::endl;
    for (size_t i = start; i < end; ++i) {
        //Current particle
        Particle& p = m_particles[i];
        
        std::vector<size_t> bonds = ::getBondPartners(bondChart, i);
        
        //Interactions
        //for (size_t j = 0; j < p.nearbyParticles.size(); ++j) {
        unsigned long int last = m_particles.size();
        if(playerParticleActivated){
            last += 1;
        }
        for (size_t j = 0; j < last; ++j) {
            bool isPlyrP = (j == m_particles.size())? true : false;
            
            //Other particle
            const Particle& q = isPlyrP? plyrP : m_particles[j];
            
            
            bool bonded = false;
            if (!isPlyrP){
                for(auto bond: bonds){
                    if(bond == j){
                        bonded = true;
                    }
                }
            }else{
                for(auto bond: plyrBondList){
                    if (bond == i){
                        bonded = true;
                    }
                }
            }
            
            //Get deltas
            float dx = q.x - p.x;
            float dy = q.y - p.y;
            if (m_wrap) {
                if (dx > m_center_x) {
                    dx -= m_width;
                } else if (dx < -m_center_x) {
                    dx += m_width;
                }
                if (dy > m_center_y) {
                    dy -= m_height;
                } else if (dy < -m_center_y) {
                    dy += m_height;
                }
            }
            float minR = m_types.MinR(p.type, q.type);
            float maxR = m_types.MaxR(p.type, q.type);
                       
            float maxR2 = maxR*maxR;
                       
           if ((std::max(dx, dy) > maxR2) && !bonded) {
               continue;
           }
            //Get distance squared
            const float r2 = dx*dx + dy*dy;
            
            
          
            
            if ((r2 > maxR2 || r2 < 0.01f) && !bonded) {
                continue;
            }
            
            if (isPlyrP){
                minR = plyrMinR;
                maxR = plyrMaxR;
            }
            
            float bondFactor = (abs(m_types.Attaract(p.type, q.type)+m_attract_mean)/(m_attract_mean+1.5*m_attract_std))/2;
            int qENum = m_types.ENumber(q.type);
            
            if(isPlyrP){
                bondFactor = (60+(rand()%40))/100.0f;
                qENum = 1;
            }
            
            int peorbsAdjust = bonded? p.e_orbitals-qENum : p.e_orbitals;
            float bondFactor2 = ((abs(peorbsAdjust+qENum)%9)/8.0f);
            
            const float bondR = 1.02*bondFactor*minR*bondFactor2; //minR/100;
            
            //            if ( rand()%100000000 < 2){
            //                std::cout << "++++++++++++++"<< std::endl;
            //                std::cout << "Electron Orbital: " << p.e_orbitals << std::endl;
            //                std::cout << "p factor: "<< abs(p.e_orbitals-8)/8.0f << std::endl;
            //                std::cout << "q factor: " <<abs(q.e_orbitals-8)/8.0f << std::endl;
            //                std::cout << "BONDS: " << bonds.size() << std::endl;
            //                std::cout << "BOND R: " << bondR << std::endl;
            //                std::cout << "MIN R: " << minR << std::endl;
            //               std::cout << "MAX R: " << maxR << std::endl;
            //                std::cout << "attract: " << m_types.Attaract(p.type, q.type)<< std::endl;
            //
            //               std::cout << "bond factor 2: " <<bondFactor2<< std::endl;
            //
            //
            //           }
            //
           
            
            //Normalize displacement
            const float r = std::sqrt(r2);
            dx /= r;
            dy /= r;
            
            bool brokeBonds = false;
            
            if (bonding && !isPlyrP){
                //Make Bond
                int p1eN = m_types.ENumber(p.type);
                int p2eN = m_types.ENumber(q.type);
                if(p.e_orbitals+p2eN <= 8 && q.e_orbitals+p1eN <= 8){
                    //override bond
                    if (!bonded &&  r < bondR){
                        
                        bonded = true;
                        bondChanges.push_back(std::make_tuple(i, j, true));
                        
                    }
                }
                
            }
            
            if (bonding && isPlyrP && plyrBonding && plyrBondList.size() < plyrBondsMax){
                if (!bonded &&  r < bondR){
                    bonded = true;
                    bondChanges.push_back(std::make_tuple(i, j, true));
                }
            }
            
            //Break Bond TODO
            if(!bonded && r<bondR/10 && rand()<1000){
                brokeBonds = true;
                for(size_t bond: bonds){
                    bondChanges.push_back(std::make_tuple(i, bond, false));
                }
                if (::isPlyrBonded(plyrBondList, i)){
                    bondChanges.push_back(std::make_tuple( m_particles.size(),i, false));
                }
                //bonds.empty();
            }
            
            if(bonded){
                
                if((r>maxR*(bondFactor2) && (isPlyrP? r>30+plyrMaxR/((((int)i*3)/m_particles.size())+1) : true)) || p.e_orbitals>8 || q.e_orbitals>8){
                    brokeBonds = true;
                    isPlyrP?  bondChanges.push_back(std::make_tuple(m_particles.size(),i, false)) :  bondChanges.push_back(std::make_tuple(i, j, false));
                    bonded = false;
                    //bondChanges.push_back(std::make_tuple(i, j, false));
//what the fuck was I doing below I'm gonna get rid of it also I should have saved a copy of this file before I started messing with it... shit.
//                    int l = 0;
//                    for(auto bond: bonds){
//                        if(bond == (int)j){
//                            //std::lock_guard<std::mutex> guard(bond_mutex);
//                            //std::cout << l << std::endl;
//                            if (l < bonds.size()){
//                                bonds.erase(bonds.begin()+l);
//                            }
//                        }
//                        l++;
//                    }
                }
            }
            
            
            //Calculate force
            
            float f = 0.0f;
            if (r > minR) {
                if (m_flat_force) {
                    f = m_types.Attaract(p.type, q.type);
                } else {
                    const float numer = 2.0f * std::abs(r - 0.5f*(maxR + minR));
                    const float denom = maxR - minR;
                    float attract = isPlyrP? plyrAttractDeNone*plyrAttractionMagnitude : m_types.Attaract(p.type, q.type);
                    f =  attract * (1.0f - numer / denom);
                }
                
            } else {
                f = R_SMOOTH*minR*(1.0f/(minR + R_SMOOTH) - 1.0f / (r + R_SMOOTH));
            }
            
            float bondPoint = 15.0f;
            
            if(isPlyrP){
                bondPoint = plyrMaxR/((((int)i*3)/m_particles.size())+1.5);
            }
            
            if (bonding && bonded){
                //if (r-bondPoint < 0){
                //f = 0.9f*(1/(1/(bondPoint+2)-1/(r+2)));//
                //f = 100*pow(bondPoint-r, 3);
                f = 0.09f*bondFactor2*(r-bondPoint);
                
            }
            
            
            /////////////////////////////////////////////////////////////
            if (!isPlyrP){
                if (dualExplosions && !bonded){
                    if (r < minR/2.5 && rand()%100>95 /*&& p.type == q.type*/){
                        //std::cout << "start type:" << (int)p.type << std::endl;
                        std::uniform_int_distribution<int>    rand_type(0, int(m_types.Size() - 1));
                        if(p.type != q.type){
                            p.type = uint8_t((int)(p.type + ((p.type+q.type)+(maxR/minR)*100))%m_types.Size());
                        }
                        //p.type = (uint8_t)rand_type(m_rand_gen);
                        
                        //std::cout << "end type:" << (int)p.type << std::endl;
                        
                    }
                }
                
                if (singleExplosions && !bonded){
                    if (r<minR/2 && rand()%100 > 98){
                        if ((int)p.type < m_types.Size()-1){
                            p.type += 1;
                        }else{
                            p.type = uint8_t(0);
                        }
                    } else if (r<minR/2.6 && rand()%1000 > 9995){
                        p.type = (uint8_t)((int)(p.type+pow(maxR/minR, 2))%m_types.Size());
                    }
                    f *= 2; //?
                }
                
                if (randomDecay){
                    if (rand()<randomDecayFactor){
                        if ((int)p.type < m_types.Size()-1){
                            try{
                                p.type = randomTypeAssignment[p.type];
                            }catch(std::exception e){
                                std::cout << "FAIL: " << p.type << std::endl;
                                std::cout << (int)p.type << std::endl;
                            }
                        }else{
                            p.type = uint8_t(0);
                        }
                    }
                }
            }
            if(brokeBonds){
                f*=3;
            }
            //Apply force
            p.vx += f * dx;
            p.vy += f * dy;
        }
    }
    return bondChanges;
}



// NEARBY PARTICLES?

//void Universe::NearbyParticlesChunk(size_t start, size_t end, float greatestV){
//    for (size_t i = start; i < end; ++i) {
//        //Current particle
//        Particle& p = m_particles[i];
//        p.nearbyParticles.clear();
//
//        //Get deltas
//        for (size_t j = 0; j < m_particles.size(); ++j) {
//            const Particle& q = m_particles[j];
//
//            float dx = q.x - p.x;
//            float dy = q.y - p.y;
//
//            if (m_wrap) {
//                if (dx > m_center_x) {
//                    dx -= m_width;
//                } else if (dx < -m_center_x) {
//                    dx += m_width;
//                }
//                if (dy > m_center_y) {
//                    dy -= m_height;
//                } else if (dy < -m_center_y) {
//                    dy += m_height;
//                }
//        }
//        //Get distance squared
//        const float r2 = dx*dx + dy*dy;
//            if (r2<greatestV*20){
//                p.nearbyParticles.push_back(q);
//            }
//            if (i == 0){
//            //std::cout << "##############" << std::endl;
//            //std::cout << p.nearbyParticles.size() << std::endl;
//            }
//            }
//    }
//}
//
//void Universe::AssignNearbyParticles(){
//    float greatestV = 0;
//    for(int i = 0; i<m_particles.size();i++){
//        float diagV = m_particles[i].vx*m_particles[i].vy;
//        if(diagV > greatestV){
//            greatestV = diagV;
//        }
//    }
//
//    size_t start = 0;
//
//    size_t end = m_particles.size();
//    size_t chunk = (numthreads>0 ? (end-start + (numthreads-1)/numthreads) : end - start);
//    std::thread tl[numthreads];
//    for( int i = 0; i<numthreads;i++){
//        size_t s = start + i*chunk;
//        size_t e = start + chunk;
//        tl[i] = std::thread{&Universe::NearbyParticlesChunk, this, s, e, greatestV};
//    }
//
//    for(int i = 0; i<numthreads;i++){
//        tl[i].join();
//    }
//}
