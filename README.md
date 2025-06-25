# GalaxySim
GalaxySim is a Python-based 3D simulation of celestial objects that showcases gravitational, thermodynamic, and luminous interactions between them.

**Features:**
1. Realistic gravitational effects using Newton's law of universal gravitation.
2. Stellar classification based on elemental composition and size.
3. Blackbody temperature and luminosity modeling.
4. Collision detection and dynamic merging of celestial bodies.
5. Interactive 3D visuals with hover-based object inspection.

**Dependencies**
Ensure the following Python modules are installed:
1. VPython – for real-time 3D visualization.
2. Mendeleev – to retrieve physical properties of elements.
3. Random – for natural variability in light placement.
4. Math – for mathematical and physical calculations.

**Installation & Usage**
1. Open your preferred Python environment (for example, IDLE, VS Code, or any IDE supporting VPython).
2. Install the required packages using the following commands:
    > pip install vpython
    > pip install mendeleev
3. Run the simulation script.
4. Enter parameters when prompted. These include the number of celestial objects and their respective name, mass, elemental composition, position, and inital velocity.
7. Once input is complete, your default web browser will launch with a real-time 3D simulation on a localhost page.
8. There, you can:
    > Pause the simulation using the on-screen button.
    > Hover over objects to view their name, classification, and temperature.
