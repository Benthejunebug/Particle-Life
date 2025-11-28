# Particle-Life
A game of life with particles: https://youtu.be/Z_zmZ23grXE

## Fork Features
This fork introduces several new mechanics and interactive elements to the simulation:

### New Mechanics
- **Player Particle**: A controllable particle that can interact with the simulation.
- **Bonding**: Particles can form bonds with each other based on proximity and type.
- **Explosions**: Trigger dual-type or single-type explosions to disrupt particle clusters.
- **Decay**: Random decay factor can be applied to particle interactions.

### Controls

#### Player Controls
- **I / J / K / L**: Move Player (Up / Left / Down / Right)
- **Num 6**: Activate/Deactivate Player Particle
- **Num 7**: Toggle Player Bonding
- **, / .**: Decrease/Increase Player Attraction
- **/**: Toggle Player Attraction Mode (None/Repel/Attract)
- **[ / ]**: Decrease/Increase Player Max Radius

#### Simulation Interaction
- **Num 2**: Toggle Dual Type Explosions
- **Num 3**: Toggle Single Type Explosions
- **Num 4**: Toggle Random Decay
- **Num 5**: Toggle Bonding
- **Left / Right**: Decrease/Increase Random Decay Factor
- **\**: Toggle Bonding Method

#### View Controls
- **Mouse Wheel**: Zoom In/Out
- **Left Click**: Track particle under cursor
- **Right Click**: Stop tracking

#### General Controls
- **Space**: Toggle Slow Motion (1 step/frame)
- **Num 0**: Fast Forward (skip drawing)
- **Num 1**: Toggle Drawing (for max speed)
- **Num 9**: Pause/Resume
- **Enter**: Randomize Particles
- **W**: Toggle Wrap-around
- **Tab**: Print Parameters
- **- / =**: Decrease/Increase Steps per Frame
