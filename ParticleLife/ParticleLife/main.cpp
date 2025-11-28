#include "Universe.h"
#include <SFML/Graphics.hpp>
#include <iostream>
#include <functional>


//
//static const int window_w = 1080;
//static const int window_h = 720;


//static const int window_w = 3200;
//static const int window_h = 1800;

//
//static const int window_w = 2000;
//static const int window_h = 2000;

//
//static const int window_w = 1920;
//static const int window_h = 1080;

static const int window_w = 1920;
static const int window_h = 1200;

static const int steps_per_frame_normal = 3;

int stepCounter = 0;
clock_t deltaTime = 0;

double clockToMilliseconds(clock_t ticks){
    // units/(units/time) => time (seconds) * 1000 = milliseconds
    return (ticks/(double)CLOCKS_PER_SEC)*1000.0;
}

int main(int argc, char *argv[]) {
    std::cout << "=========================================================" << std::endl;
    std::cout << std::endl;
    std::cout << "               Welcome to Particle Life" << std::endl;
    std::cout << std::endl;
    std::cout << "  This is a particle-based game of life simulation based" << std::endl;
    std::cout << "on random attraction and repulsion between all particle" << std::endl;
    std::cout << "classes.  For more details about how this works and other" << std::endl;
    std::cout << "fun projects, check out my YouTube channel 'CodeParade'." << std::endl;
    std::cout << std::endl;
    std::cout << "=========================================================" << std::endl;
    std::cout << std::endl;
    std::cout << "Controls:" << std::endl;
    std::cout << "         'B' - Randomize (Balanced)" << std::endl;
    std::cout << "         'C' - Randomize (Chaos)" << std::endl;
    std::cout << "         'D' - Randomize (Diversity)" << std::endl;
    std::cout << "         'F' - Randomize (Frictionless)" << std::endl;
    std::cout << "         'G' - Randomize (Gliders)" << std::endl;
    std::cout << "         'H' - Randomize (Homogeneity)" << std::endl;
    std::cout << "         'L' - Randomize (Large Clusters)" << std::endl;
    std::cout << "         'M' - Randomize (Medium Clusters)" << std::endl;
    std::cout << "         'Q' - Randomize (Quiescence)" << std::endl;
    std::cout << "         'S' - Randomize (Small Clusters)" << std::endl;
    std::cout << "         'W' - Toggle Wrap-Around" << std::endl;
    std::cout << "       Enter - Keep rules, but re-seed particles" << std::endl;
    std::cout << "       Space - Hold for slow-motion" << std::endl;
    std::cout << "         Tab - Print current parameters to console" << std::endl;
    std::cout << "  Left Click - Click a particle to follow it" << std::endl;
    std::cout << " Right Click - Click anywhere to unfollow particle" << std::endl;
    std::cout << "Scroll Wheel - Zoom in/out" << std::endl;
    std::cout << std::endl;
    //system("pause");
    
    //Create the universe of particles
    Universe universe(9, 400, window_w, window_h);
    universe.ReSeed(-0.02f, 0.06f, 0.0f, 20.0f, 20.0f, 70.0f, 0.05f, false);
    
    //set to my test
    universe.SetPopulation(10, 2000);
    universe.ReSeed(0.0065f, 0.09f, 0.0f, 20.0f, 40.0f, 110.0f, 0.093134f, false);
    //universe.ToggleWrap();
    
    //Camera settings
    float cam_x = float(window_w/2);
    float cam_y = float(window_h/2);
    float cam_zoom = 1.0f;
    float cam_x_dest = cam_x;
    float cam_y_dest = cam_y;
    float cam_zoom_dest = cam_zoom;
    int32_t last_scroll_time = 0;
    int track_index = -1;
    int steps_per_frame = steps_per_frame_normal;
    
    //GL settings
    sf::ContextSettings settings;
    settings.depthBits = 24;
    settings.stencilBits = 8;
    settings.antialiasingLevel = 5;
    //settings.majorVersion = 3;
    //settings.minorVersion = 0;
    
    //Create the window
    const sf::VideoMode screenSize = sf::VideoMode(window_w, window_h, 24);
    sf::RenderWindow window(screenSize, "Particles", sf::Style::Resize | sf::Style::Close, settings);
    window.setFramerateLimit(30);
    window.setVerticalSyncEnabled(true);
    window.setActive(false);
    window.requestFocus();
    sf::Clock clock;
    
    //set view to whole universe in case the window wouldn't normally fit on the screen had to be added when moving to
    // sfml 2.6.1 as views may not have been a thing before otherwise he would have used that to do the zoom and such
    // but this works so whatever
    sf::FloatRect visibleArea(0, 0, window_w, window_h);
    window.setView(sf::View(visibleArea));

    
    //Ben's Settings
    bool skipDraw = false;
    bool dontDraw = false;
    bool pause = false;
    
    //Main Loop
    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window.close();
                break;
            } else if (event.type == sf::Event::Resized){
                sf::FloatRect visibleArea(0, 0, window_w, window_h);
                window.setView(sf::View(visibleArea));
            } else if (event.type == sf::Event::KeyPressed) {
                const sf::Keyboard::Key keycode = event.key.code;
                if (keycode == sf::Keyboard::Escape) {
                    window.close();
                    break;
                } else if (keycode == sf::Keyboard::P) { //Lots Medium Clusters
                    universe.SetPopulation(10, 1000);
                    universe.ReSeed(0.045f, 0.03f, 0.0f, 30.0f, 30.0f, 80.0f, 0.03f, false);
                } else if (keycode == sf::Keyboard::T) { //Ben's Test
                    universe.SetPopulation(10, 3000);
                    universe.ReSeed(0.0065f, 0.09f, 0.0f, 20.0f, 40.0f, 70.0f, 0.1099134f, false);
                } else if (keycode == sf::Keyboard::Y) { //Ben's Test 2 ---- REALLY COOL but CONFUSING
                    universe.SetPopulation(100, 2500);
                    universe.ReSeed(0.025f, 0.832666f, 0.0f, 10.0f, 40.0f, 520.0f, 0.8043134f, false);
                } else if (keycode == sf::Keyboard::U) { //Ben's Test 2 ---- REALLY COOL but CONFUSING --->>> current
                    universe.SetPopulation(6, 3000);
                    universe.ReSeed(0.025f, 0.082666f, 0.0f, 10.0f, 40.0f, 1000.0f, 0.9843134f, false);
                } else if (keycode == sf::Keyboard::I) { //Ben's Test 3 ---- REALLY COOL but CONFUSING --->>> current
                    universe.SetPopulation(10, 2000);
                    universe.ReSeed(0.027f, 0.12666f, 0.0f, 10.0f, 40.0f, 2000.0f, 0.7543134f, false);
                } else if (keycode == sf::Keyboard::B) { //Balanced
                    universe.SetPopulation(9, 2000);
                    universe.ReSeed(-0.02f, 0.08f, 0.0f, 20.0f, 20.0f, 70.0f, 0.05f, false);
                } else if (keycode == sf::Keyboard::C) { //Chaos
                    universe.SetPopulation(6, 2000);
                    universe.ReSeed(0.02f, 0.04f, 0.0f, 30.0f, 30.0f, 100.0f, 0.01f, false);
                } else if (keycode == sf::Keyboard::D) { //Diversity
                    universe.SetPopulation(16, 3000);
                    universe.ReSeed(0.01f, 0.06f, 0.0f, 20.0f, 20.0f, 60.0f, 0.05f, false);
                } else if (keycode == sf::Keyboard::F) { //Frictionless
                    universe.SetPopulation(4, 500);
                    universe.ReSeed(0.01f, 0.005f, 10.0f, 10.0f, 10.0f, 60.0f, 0.0f, true);
                } else if (keycode == sf::Keyboard::G) { //Gliders
                    universe.SetPopulation(6, 1000);
                    universe.ReSeed(0.0f, 0.06f, 0.0f, 20.0f, 10.0f, 50.0f, 0.1f, true);
                } else if (keycode == sf::Keyboard::H) { //Homogeneity
                    universe.SetPopulation(4, 1000);
                    universe.ReSeed(0.0f, 0.1f, 10.0f, 10.0f, 10.0f, 80.0f, 0.05f, true);
                } else if (keycode == sf::Keyboard::N) { //Large Clusters
                    universe.SetPopulation(8, 2000);
                    universe.ReSeed(0.025f, 0.04f, 0.0f, 30.0f, 30.0f, 100.0f, 0.3f, false);
                } else if (keycode == sf::Keyboard::M) { //Medium Clusters
                    universe.SetPopulation(8, 1000);
                    universe.ReSeed(0.02f, 0.05f, 0.0f, 20.0f, 20.0f, 50.0f, 0.05f, false);
                } else if (keycode == sf::Keyboard::Q) { //Quiescence
                    universe.SetPopulation(6, 300);
                    universe.ReSeed(-0.02f, 0.04f, 10.0f, 20.0f, 20.0f, 60.0f, 0.2f, false);
                } else if (keycode == sf::Keyboard::S) { //Small Clusters
                    universe.SetPopulation(6, 1000);
                    universe.ReSeed(-0.005f, 0.01f, 10.0f, 10.0f, 20.0f, 50.0f, 0.01f, false);
                } else if (keycode == sf::Keyboard::W) {
                    universe.ToggleWrap();
                } else if (keycode == sf::Keyboard::Enter) {
                    universe.SetRandomParticles();
                } else if (keycode == sf::Keyboard::Tab) {
                    universe.PrintParams();
                } else if (keycode == sf::Keyboard::Space) {
                    steps_per_frame = 1;
                    skipDraw = false;
                } else if (keycode == sf::Keyboard::Equal) {
                    steps_per_frame += 10;
                } else if (keycode == sf::Keyboard::Num0) {
                    steps_per_frame += 100;
                    skipDraw = true;
                } else if (keycode == sf::Keyboard::Hyphen) {
                    steps_per_frame -= 10;
                } else if (keycode == sf::Keyboard::Num1) {
                    dontDraw = !dontDraw;
                    if(dontDraw){
                        steps_per_frame = 100;
                    }else{
                        steps_per_frame = steps_per_frame_normal;
                    }
                } else if (keycode == sf::Keyboard::Num2) {
                    universe.ToggleDualExplosions();
                } else if (keycode == sf::Keyboard::Num3){
                    universe.TogglSingleExplosions();
                } else if (keycode == sf::Keyboard::Num4){
                    universe.TogglRandomDecay();
                } else if (keycode == sf::Keyboard::Num5){
                    universe.TogglBonding();
                } else if (keycode == sf::Keyboard::Num6){
                    universe.ActivatePlayerParticle();
                } else if (keycode == sf::Keyboard::Num7){
                    universe.TogglPlyrBonding();
                } else if (keycode == sf::Keyboard::Num9){
                    pause = !pause;
                } else if (keycode == sf::Keyboard::Left){
                    universe.DecreaseRandomDecayFactor();
                }   else if (keycode == sf::Keyboard::Right){
                    universe.IncreaseRandomDecayFactor();
                }else if (keycode == sf::Keyboard::Comma){
                    universe.ChangePlayerAttract(false);
                }else if (keycode == sf::Keyboard::Period){
                    universe.ChangePlayerAttract(true);
                }else if (keycode == sf::Keyboard::Slash){
                    universe.TrogglePlayerAttract();
                }else if (keycode == sf::Keyboard::LBracket){
                    universe.ChangePlyrMaxR(false);
                }else if (keycode == sf::Keyboard::RBracket){
                    universe.ChangePlyrMaxR(true);
                }
                
                if (keycode == sf::Keyboard::I){
                    universe.ChangePlayerDirection(universe.PLYRUP, true);
                }
                if (keycode == sf::Keyboard::J){
                    universe.ChangePlayerDirection(universe.PLYRLF, true);
                }
                if (keycode == sf::Keyboard::K){
                    universe.ChangePlayerDirection(universe.PLYRDN, true);
                }
                if (keycode == sf::Keyboard::L){
                    universe.ChangePlayerDirection(universe.PLYRRT, true);
                }
                
            } else if (event.type == sf::Event::KeyReleased) {
                
                const sf::Keyboard::Key keycode = event.key.code;
                
                if (keycode == sf::Keyboard::Space) {
                    steps_per_frame = steps_per_frame_normal;
                } else if (keycode == sf::Keyboard::I){
                    universe.ChangePlayerDirection(universe.PLYRUP, false);
                } else if (keycode == sf::Keyboard::J){
                    universe.ChangePlayerDirection(universe.PLYRLF, false);
                }else if (keycode == sf::Keyboard::K){
                    universe.ChangePlayerDirection(universe.PLYRDN, false);
                }else if (keycode == sf::Keyboard::L){
                    universe.ChangePlayerDirection(universe.PLYRRT, false);
                }else if (keycode == sf::Keyboard::BackSlash){
                    universe.TroggleNewBondingCode();
                }
                
            } else if (event.type == sf::Event::MouseWheelMoved) {
                cam_zoom_dest *= std::pow(1.1f, event.mouseWheel.delta);
                cam_zoom_dest = std::max(std::min(cam_zoom_dest, 10.0f), 1.0f);
                const int32_t cur_time = clock.getElapsedTime().asMilliseconds();
                if (cur_time - last_scroll_time > 300) {
                    //Only update position if scroll just started
                    const sf::Vector2i mouse_pos = sf::Mouse::getPosition(window);
                    universe.ToCenter(mouse_pos.x, mouse_pos.y, cam_x_dest, cam_y_dest);
                }
                last_scroll_time = cur_time;
            } else if (event.type == sf::Event::MouseButtonPressed) {
                if (event.mouseButton.button == sf::Mouse::Left) {
                    track_index = universe.GetIndex(event.mouseButton.x, event.mouseButton.y);
                } else if (event.mouseButton.button == sf::Mouse::Right) {
                    track_index = -1;
                }
            }
        }
        
        //Apply zoom
        if (track_index >= 0) {
            cam_x_dest = universe.GetParticleX(track_index);
            cam_y_dest = universe.GetParticleY(track_index);
        }
        cam_x = cam_x*0.9f + cam_x_dest*0.1f;
        cam_y = cam_y*0.9f + cam_y_dest*0.1f;
        cam_zoom = cam_zoom*0.8f + cam_zoom_dest*0.2f;
        universe.Zoom(cam_x, cam_y, cam_zoom);
        
        if (pause){
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
            continue;
        }
        //Apply physics and draw
        //window.clear(sf::Color::White);
        window.clear();
        for (int i = 0; i < steps_per_frame; ++i) {
            const float opacity = pow(float(i + 1) / float(steps_per_frame), 2);
            clock_t beginFrame = std::clock();
            universe.Step();
            clock_t endFrame = std::clock();
            
            deltaTime += endFrame - beginFrame;
            stepCounter ++;
            
            if( clockToMilliseconds(deltaTime)>1000.0){
                std::cout << "SPS: " << stepCounter << std::endl;
                stepCounter = 0;
                deltaTime = 0;
            }
            if(!skipDraw && !dontDraw){
                universe.Draw(window, opacity);
            }
        }
        if(skipDraw && !dontDraw){
            universe.Draw(window, 1);
        }
        //Flip the screen buffer
        window.display();
    }
    
    return 0;
}
