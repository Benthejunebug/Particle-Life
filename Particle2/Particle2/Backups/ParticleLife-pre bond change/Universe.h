#pragma once
#include "Particles.h"
#include <random>
#include <iostream>
#include <vector>
#include <thread>
#include <algorithm>
#include <mutex>


class Universe {
public:
    
    const int PLYRNONE = 0;
    const int PLYRUP = 2;
    const int PLYRDN = 3;
    const int PLYRLF = 5;
    const int PLYRRT = 7;
    
  Universe(size_t num_types, size_t num_particles, int width, int height);

  void ReSeed(float attract_mean, float attract_std, float minr_lower, float minr_upper,
    float maxr_lower, float maxr_upper, float friction, bool flat_force);

  void SetPopulation(size_t num_types, size_t num_particles);
  void SetSize(float width, float height) { m_width = width; m_height = height; }
  void SetRandomTypes();
  void SetRandomParticles();
  void SetRandomTypeAllocation();
  void ToggleWrap() { m_wrap = !m_wrap; }
    void ToggleDualExplosions(){dualExplosions = !dualExplosions;
        std::cout << "#############" << std::endl;
        if (dualExplosions){
            std::cout << "DUAL TYPE EXPLOSIONS ON" << std::endl;
        }else{
            std::cout << "DUAL TYPE EXPLOSIONS OFF" << std::endl;
        }
    }
    void TogglSingleExplosions(){singleExplosions = !singleExplosions;
        std::cout << "#############" << std::endl;
        if (singleExplosions){
            std::cout << "SINGLE TYPE EXPLOSIONS ON" << std::endl;
        }else{
            std::cout << "SINGLE TYPE Explosions OFF" << std::endl;
        }
    }
    void TogglRandomDecay(){randomDecay = !randomDecay;
        std::cout << "#############" << std::endl;
        if (randomDecay){
            std::cout << "RANDOM DECAY ON" << std::endl;
        }else{
            std::cout << "RANDOM DECAY OFF" << std::endl;
        }
    }
    void IncreaseRandomDecayFactor(){
        randomDecayFactor += 100;
        std::cout << "#############" << std::endl;
        std::cout << "RANDOM DECAY FACTOR (+5): " << randomDecayFactor << std::endl;
    }
    void DecreaseRandomDecayFactor(){
        randomDecayFactor -= 100;
        std::cout << "#############" << std::endl;
        std::cout << "RANDOM DECAY FACTOR (-5): " << randomDecayFactor << std::endl;
    }
    void TogglBonding(){bonding = !bonding;
        std::cout << "#############" << std::endl;
        if (bonding){
            std::cout << "BONDING ON" << std::endl;
        }else{
            std::cout << "BONDING OFF" << std::endl;
        }
    }

    void ActivatePlayerParticle(){playerParticleActivated = !playerParticleActivated;
        plyrBondList.clear();
        std::cout << "#############" << std::endl;
        if (playerParticleActivated){
            std::cout << "PLAYER PARTICLE ACTIVATED" << std::endl;
        }else{
            std::cout << "PLAYER PARTICLE DEACTIVATED" << std::endl;
        }
    }
    
    bool SameDirection(int directionCode, int directionToCheck){
        if(directionCode == 1){return false;}
        return (directionCode%directionToCheck == 0);
    }
    
    void ChangePlyrMaxR(bool increaseDec){
        plyrMaxR += (increaseDec? 5 : -5);
        plyrMaxR = abs(plyrMaxR);
        std::cout << "#############" << std::endl;
        std::cout << "PLAYER MAX-R " << (increaseDec? "INCREASE: " : "DECREASE: ") << plyrMaxR << std::endl;
    }
    
    void ChangePlayerAttract(bool increaseDec){
        plyrAttractionMagnitude += (increaseDec? 1 : -1)*plyrAttractionIncrement;
        plyrAttractionMagnitude = abs(plyrAttractionMagnitude);
        std::cout << "#############" << std::endl;
        std::cout << "PLAYER ATTRACTION MAGNITUDE " << (increaseDec? "INCREASE: " : "DECREASE: ") << plyrAttractionMagnitude << std::endl;
        //SetPlyrMaxR();
    }
    
    
    void ClearPlayerDirection(){
        plyrDirection = 1;
    }
    
    void ChangePlayerDirection(int d, bool add){
        if(((float)plyrDirection/(float)d) >= 1.0f && add){
            plyrDirection = plyrDirection%d;
        }
        plyrDirection = (!add? plyrDirection/d : plyrDirection*d);
        if (plyrDirection > 35 || plyrDirection == 0){
            plyrDirection = 1;
        }
    }
    
    //threeway toggle
    void TrogglePlayerAttract(){
        plyrAttractDeNone = (plyrAttractDeNone+2)%3-1;
        std::string message;
        if (plyrAttractDeNone == 0){
            message = "NONE";
        }else if(plyrAttractDeNone == -1){
            message = "REPEL";
        }else{
            message = "ATTRACT";
        }
        std::cout << "PLAYER ATTRACT?: " << message << std::endl;
    }
    
    void TogglPlyrBonding(){plyrBonding = !plyrBonding;
        std::cout << "#############" << std::endl;
        if (plyrBonding){
            std::cout << "PLAYER BONDING ON" << std::endl;
        }else{
            plyrBondList.clear();
            std::cout << "PLAYER BONDING OFF" << std::endl;
        }
    }
    
//    void SetPlyrMaxR(){
//        plyrMaxR = ((m_center_x - plyrMaxR_lower*0.6)*plyrAttractionMagnitude/5)+plyrMaxR_lower*0.8;
//    }
    
  void Step();
    void AssignNearbyParticles();
    void NearbyParticlesChunk(size_t start, size_t end, float greatestV);
    void chunkInteractions(size_t start, size_t end);
    void chunkUpdate(size_t start, size_t end);
    void updateParticle(Particle& p);
  void Draw(sf::RenderWindow& window, float opacity) const;
  void Zoom(float cx, float cy, float zoom);

  int GetIndex(int x, int y) const;
  float GetParticleX(int index) const;
  float GetParticleY(int index) const;
  void ToCenter(int x, int y, float& cx, float& cy) const;

  void PrintParams() const;
    void PrintTypeCount();

private:
    float m_center_x;
    float m_center_y;
    float m_zoom;
    bool m_wrap;
    int numthreads;
    std::vector<std::thread> th;
    int nearbyCounter;
    bool dualExplosions;
    bool singleExplosions;
    bool randomDecay;
    int randomDecayFactor;
    std::vector<int> randomTypeAssignment;
    std::vector<int> RandomizeTypes(int size);
    bool CheckTypeAssignment(std::vector<int>);
    
    bool bonding;
    std::vector<std::vector<int>> bondChart;
    std::mutex bond_mutex;
    
    
    bool bonded(int p);
    void addBond(int p1I, int p2I);
    void removeBond(int p1I, int p2I);
    std::vector<int> getBondPartners(int particleIndex);
    int getBondPartnersCount(int particleIndex);
    
    // Player Particle implementation
    bool playerParticleActivated = false;
    sf::Color plyrColor = sf::Color::White;
    int plyrAttractDeNone = 0;
    Particle plyrP;
    float plyrRadius;
    float plyrMinR;
    float plyrMaxR;
    //float plyrMaxR_lower;
    float plyrAttractionMagnitude;
    float plyrAttractionIncrement;
    float plyrV;
    int plyrDirection;
    bool plyrBonding;
    std::vector<int> plyrBondList = {};
    void addPlyrBond(int p);
    void removePlyrBond(int p);
    bool isPlyrBonded(int p);
    int plyrBondsMax;
    
    std::vector<Particle> m_particles;
    ParticleTypes m_types;
    float m_width;
    float m_height;

  std::mt19937 m_rand_gen;

  //Random settings
  float m_attract_mean;
  float m_attract_std;
  float m_minr_lower;
  float m_minr_upper;
  float m_maxr_lower;
  float m_maxr_upper;
  float m_friction;
  bool m_flat_force;
};
