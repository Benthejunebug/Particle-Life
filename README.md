# Particle-Life Fork

A fork of [CodeParade's Particle-Life](https://github.com/HackerPoet/Particle-Life) with enhanced interactive features.

## What's Different From Upstream

This fork introduces several experimental mechanics and interactive elements:

### New Mechanics

- **Player Particle**: A controllable white particle that can interact with the simulation
- **Bonding System**: Particles can form bonds with each other based on proximity and attraction rules
  - Three bonding methods: Bond Chart, Particle Bond List, and Bond Table
  - Bonds persist and affect particle movement
- **Explosion Effects**: 
  - Dual-type explosions that disrupt clusters of two particle types
  - Single-type explosions affecting one particle type
- **Random Decay**: Configurable decay factor that adds entropy to the simulation

### Controls

#### Player Particle Controls
- **I / J / K / L**: Move Player Particle (Up / Left / Down / Right)
- **6**: Toggle Player Particle on/off
- **7**: Toggle Player Bonding
- **, / .**: Decrease/Increase Player Attraction Magnitude
- **/**: Cycle Player Attraction Mode (None → Repel → Attract)
- **[ / ]**: Decrease/Increase Player Max Radius

#### Simulation Controls
- **2**: Toggle Dual Type Explosions
- **3**: Toggle Single Type Explosions
- **4**: Toggle Random Decay
- **5**: Toggle Bonding System
- **Left / Right Arrow**: Decrease/Increase Random Decay Factor
- **\\**: Cycle Bonding Method

#### View & Performance
- **Mouse Wheel**: Zoom In/Out (focused on cursor position)
- **Left Click**: Track particle under cursor
- **Right Click**: Stop tracking
- **Space**: Toggle Slow Motion (hold for 1 step/frame, release for normal speed)
- **0**: Fast Forward mode (skips drawing for max speed)
- **1**: Toggle drawing entirely (for maximum simulation speed)
- **9**: Pause/Resume simulation
- **- / =**: Decrease/Increase steps per frame

#### Presets & Utilities
- **G**: Gliders preset
- **H**: Homogeneity preset
- **N**: Large Clusters preset
- **M**: Medium Clusters preset
- **Q**: Quiescence preset
- **S**: Small Clusters preset
- **W**: Toggle wrap-around boundary
- **Enter**: Randomize particles
- **Tab**: Print current parameters to console

## Technical Notes

- Multi-threaded particle interactions for improved performance
- Three bonding implementation methods for experimentation
- Real-time SPS (Steps Per Second) counter displayed in terminal
- Configurable player particle parameters (attraction, radius, bonding behavior)

## Building

This fork requires SFML 2.5+ frameworks.

### macOS (CLI / CMake)

1.  Ensure SFML is installed in `/Library/Frameworks`.
2.  If you encounter "damaged" framework errors, remove the quarantine attributes:
    ```bash
    xattr -r -d com.apple.quarantine /Library/Frameworks/sfml-*.framework
    ```
3.  Build using CMake:
    ```bash
    cmake .
    make
    ./ParticleLife
    ```

### Xcode

The project includes Xcode project files for macOS in the `particle.xcodeproj` directory.

---

**Upstream**: [HackerPoet/Particle-Life](https://github.com/HackerPoet/Particle-Life)  
**Video**: https://youtu.be/Z_zmZ23grXE
