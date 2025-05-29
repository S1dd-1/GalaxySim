from math import *
from vpython import *
from mendeleev import *
from random import *

scene = canvas(title="Galaxy Simulation")
scene.autoscale = True
scene.lights = []

G = 6.674e-11  # Gravitational constant
c = 299792458 # Speed of light
hbar = 1.055e-34 # Reduced plank constant
sigma = 5.67e-8 # Stefan-Boltzmann constant
kb = 1.38e-23 # Boltzmann constant
Wi = 2.897e-3 # Wien's constant
Msun = 1.989e30 # Mass of Sun

abc = {}  # Main Dictionary

def Offset(R):
    return (uniform(-0.5,0.5))*R

def RGB(W):
    r = e**(-0.5*((W-610)/30)**2)
    g = e**(-0.5*((W-545)/25)**2)
    b = e**(-0.5*((W-450)/20)**2)
    return vector(r, g, b)

def Normalize(x,vecX,y,vecY,z,vecZ):
    total = x+y+z
    if total == 0:
        return vector(0.2,0.2,0.2)
    x /= total
    y /= total
    z /= total
    return vector(
        x * vecX.x + y * vecY.x + z * vecZ.x,
        x * vecX.y + y * vecY.y + z * vecZ.y,
        x * vecX.z + y * vecY.z + z * vecZ.z
    )

def Motionise(n,p,dt):
    ap = abc[n]["Obj_Fnet"][f"Fnet{p}"] / abc[n]["Obj_Mass"]
    abc[n]["Obj_Vel"][f"v{p}"] += ap * dt
    abc[n]["Obj_Position"][f"{p}"] += abc[n]["Obj_Vel"][f"v{p}"] * dt

def EquiTemp(comp, flux):
    A = 0.3 + 0.2 * (comp.get('Fe', 0)) - 0.1 * (comp.get('C', 0))
    A = max(0.0, min(A, 1.0))

    G = 1.0 + 0.3 * (comp.get('C', 0))
    G = max(1.0, G)

    if flux == 0:
        return 2.7

    T = G * (((1 - A) * flux) / (4 * sigma)) ** 0.25
    return T

def Densiter(ele):
    netEleDense = 0
    netPerc = 0
    for j in ele:
        eleDense = element(j).density * 1000
        netEleDense += eleDense * ele[j]
        netPerc += ele[j]

    netDensity = netEleDense / netPerc
    return netDensity

def distance(i,j):
    rx = abc[j]["Obj_Position"]["x"] - abc[i]["Obj_Position"]["x"]
    ry = abc[j]["Obj_Position"]["y"] - abc[i]["Obj_Position"]["y"]
    rz = abc[j]["Obj_Position"]["z"] - abc[i]["Obj_Position"]["z"]
    r = sqrt(rx**2 + ry**2 + rz**2)

    return rx, ry, rz, r

def classify(i):
    mass = abc[i]["Obj_Mass"]
    radi = abc[i]["Obj_Radius"][0]

    if radi <= abc[i]["Obj_Radius"][1]:
        abc[i]["Obj_Type"] = "BlackHole"
        abc[i]["Obj_Temp"] = (hbar*(c**3))/(8*pi*G*mass*kb)
        abc[i]["Obj_Color"] = vector(0.1,0.1,0.1)

    elif abc[i]["Obj_Comp"].get("H", 0) >= 0.4 and abc[i]["Obj_Comp"].get("He", 0) > 0.1:
        if abc[i]["Obj_Mass"] > 0.008*Msun:
            abc[i]["Obj_Type"] = "Star"
            abc[i]["Obj_Temp"] = 5772 * (mass / Msun)**0.54
            abc[i]["Obj_Lumin"]= 4*pi*(radi**2)*sigma*(abc[i]["Obj_Temp"]**4)
            lightMagnitude = int(round(0.0748*(abc[i]["Obj_Temp"]**0.445)))
            waveL = max(400,min(Wi/abc[i]["Obj_Temp"],700))
            abc[i]["Obj_Color"]= RGB(waveL)

            starlight = []
            for j in range(lightMagnitude):
                starlight.append(local_light(pos=vector(abc[i]["Obj_Position"]["x"] + Offset(radi), abc[i]["Obj_Position"]["y"] + Offset(radi), abc[i]["Obj_Position"]["z"] + Offset(radi)), color=abc[i]["Obj_Color"]))
            abc[i]["Obj_Lights"] = starlight

            abc[i]["Obj_Sphere"] = sphere(pos=vector(abc[i]["Obj_Position"]["x"], abc[i]["Obj_Position"]["y"], abc[i]["Obj_Position"]["z"]), radius=radi, color=abc[i]["Obj_Color"], make_trail=True, retain=1000, emissive = True)
            
        elif abc[i]["Obj_Comp"].get("O", 0) > 0.15:
            abc[i]["Obj_Type"] = "IceGiant"
            abc[i]["Obj_Color"] = vector(0.2, 0.5, 1.0)
            abc[i]["Obj_Sphere"] = sphere(pos=vector(abc[i]["Obj_Position"]["x"], abc[i]["Obj_Position"]["y"], abc[i]["Obj_Position"]["z"]), radius=radi, color=abc[i]["Obj_Color"], make_trail=True, retain=1000, emissive = False)
                        
        else:
            abc[i]["Obj_Type"] = "GasGiant"
            abc[i]["Obj_Color"] = vector(1.0, 0.9, 0.6)
            abc[i]["Obj_Color"] = Normalize(abc[i]["Obj_Comp"].get("C", 0), vector(0.2, 0.5, 1.0), abc[i]["Obj_Comp"].get("N", 0), vector(1.0, 1.0, 0.7), abc[i]["Obj_Comp"].get("O", 0),vector(0.5, 0.7, 1.0))
            abc[i]["Obj_Sphere"] = sphere(pos=vector(abc[i]["Obj_Position"]["x"], abc[i]["Obj_Position"]["y"], abc[i]["Obj_Position"]["z"]), radius=radi, color=abc[i]["Obj_Color"], make_trail=True, retain=1000, emissive = False)

    else:
        if radi <= 2370e3:
            abc[i]["Obj_Type"] = "DwarfPlanet"
            abc[i]["Obj_Color"] = vector(0.8, 0.8, 0.9)
            abc[i]["Obj_Sphere"] = sphere(pos=vector(abc[i]["Obj_Position"]["x"], abc[i]["Obj_Position"]["y"], abc[i]["Obj_Position"]["z"]), radius=radi, color=abc[i]["Obj_Color"], make_trail=True, retain=1000, emissive = False)
        else:
            abc[i]["Obj_Type"] = "Planet"
            abc[i]["Obj_Color"] = Normalize(abc[i]["Obj_Comp"].get("O", 0), vector(0.4, 0.7, 0.9), abc[i]["Obj_Comp"].get("Fe", 0), vector(0.8, 0.4, 0.3), abc[i]["Obj_Comp"].get("C", 0),vector(1.0, 0.8, 0.4))
            abc[i]["Obj_Sphere"] = sphere(pos=vector(abc[i]["Obj_Position"]["x"], abc[i]["Obj_Position"]["y"], abc[i]["Obj_Position"]["z"]), radius=radi, color=abc[i]["Obj_Color"], make_trail=True, retain=1000, emissive = False)

ObjectNum = int(input("Enter number of objects: "))

for i in range(ObjectNum):
    name = input("Enter the celestial object's name: ")
    mass = float(input(f"Enter the mass of {name}: "))
    
    elenum = int(input(f"Enter the number of elements the {name} is comprised of: "))
    netEleDense = 0
    netPerc = 0
    ele = {}
    for j in range(elenum):
        eleSym = input("Enter the element symbol: ").title()
        elePerc = float(input("Enter the percentage composition: "))/100
        ele[eleSym] = elePerc

    radi = ((3 * mass) / (4 * pi * Densiter(ele))) ** (1 / 3)
    schwarzRadi = (2*G*mass)/(c**2)
    x = float(input(f"Enter x-coordinate of {name}: "))
    y = float(input(f"Enter y-coordinate of {name}: "))
    z = float(input(f"Enter z-coordinate of {name}: "))
    vx = float(input(f"Enter x-component of {name}'s velocity: "))
    vy = float(input(f"Enter y-component of {name}'s velocity: "))
    vz = float(input(f"Enter z-component of {name}'s velocity: "))

    abc[i] = {
        "Obj_Name": name,
        "Obj_Comp": ele,
        "Obj_Mass": mass,
        "Obj_Density": Densiter(ele),
        "Obj_Radius": [radi,schwarzRadi],
        "Obj_Position": {"x": x, "y": y, "z": z},
        "Obj_Vel": {"vx": vx, "vy": vy, "vz": vz},
        "Obj_Fnet": {"Fnetx": 0, "Fnety": 0, "Fnetz": 0},
        "Obj_Lights" : []
    }

    classify(i)
    
def Pause():
    global running  
    running = not running  
    if running == True:
        pausebut.text = "Pause"
    elif running == False:
        pausebut.text = "Play"

pausebut = button(bind=Pause, text="Pause")

running = True
dt = 0.01
abcKeys = list(abc.keys())

while True:
    rate(60)
    if running:
        ToRemove = []
        for i in abcKeys:
            if i in ToRemove:
                continue
            else:
                radi = abc[i]["Obj_Radius"][0]
                mass1 = abc[i]["Obj_Mass"]
                for j in abc:
                    if j in ToRemove:
                        continue
                    else:
                        if j > i:
                            mass2 = abc[j]["Obj_Mass"]
                            rx, ry, rz, r = distance(i, j)

                            if r <= radi+abc[j]["Obj_Radius"][0]:
                                #Collision
                                abc[i]["Obj_Name"] = abc[i]["Obj_Name"] + "-" + abc[j]["Obj_Name"]
                                abc[i]["Obj_Mass"] = mass1 + mass2
                                abc[i]["Obj_Vel"]["vx"] = (mass1*abc[i]["Obj_Vel"]["vx"] + mass2*abc[j]["Obj_Vel"]["vx"])/(mass1+mass2)
                                abc[i]["Obj_Vel"]["vy"] = (mass1*abc[i]["Obj_Vel"]["vy"] + mass2*abc[j]["Obj_Vel"]["vy"])/(mass1+mass2)
                                abc[i]["Obj_Vel"]["vz"] = (mass1*abc[i]["Obj_Vel"]["vz"] + mass2*abc[j]["Obj_Vel"]["vz"])/(mass1+mass2)

                                for k in {**abc[i]["Obj_Comp"],**abc[j]["Obj_Comp"]}:
                                    abc[i]["Obj_Comp"][k] = (mass1*abc[i]["Obj_Comp"].get(k,0) + mass2*abc[j]["Obj_Comp"].get(k,0))/(mass1+mass2)
                                abc[i]["Obj_Density"] = Densiter(abc[i]["Obj_Comp"])
                                abc[i]["Obj_Radius"][0] = ((3*abc[i]["Obj_Mass"]) / (4*pi*abc[i]["Obj_Density"])) ** (1 / 3)
                                abc[i]["Obj_Radius"][1] = (2*G*abc[i]["Obj_Mass"])/(c**2)

                                classify(i)
                                ToRemove.append(j)
                            else:
                                F = (G * mass2 * mass1) / (r**2)
                                Fx = F * (rx / r)
                                Fy = F * (ry / r)
                                Fz = F * (rz / r)

                                abc[i]["Obj_Fnet"]["Fnetx"] += Fx
                                abc[j]["Obj_Fnet"]["Fnetx"] -= Fx
                                abc[i]["Obj_Fnet"]["Fnety"] += Fy
                                abc[j]["Obj_Fnet"]["Fnety"] -= Fy
                                abc[i]["Obj_Fnet"]["Fnetz"] += Fz
                                abc[j]["Obj_Fnet"]["Fnetz"] -= Fz

            for i in abcKeys:
                if i in ToRemove:
                    continue
                else:
                    radi = abc[i]["Obj_Radius"][0]

                    Motionise(i,"x",dt)
                    Motionise(i,"y",dt)
                    Motionise(i,"z",dt)

                    abc[i]["Obj_Sphere"].pos = vector(abc[i]["Obj_Position"]["x"], abc[i]["Obj_Position"]["y"], abc[i]["Obj_Position"]["z"])
                    for j in abc[i]["Obj_Lights"]:
                        j.pos = vector(abc[i]["Obj_Position"]["x"] + Offset(radi), abc[i]["Obj_Position"]["y"]+ Offset(radi), abc[i]["Obj_Position"]["z"] + Offset(radi))

            for i in abc:
                if i in ToRemove:
                    continue
                else:
                    if abc[i]["Obj_Type"] in ["Planet", "DwarfPlanet", "IceGiant", "GasGiant"]:
                        netflux = 0
                        for j in abc:
                            if abc[j]["Obj_Type"] == "Star":
                                rx, ry, rz, r = distance(i, j)
                                flux = abc[j]["Obj_Lumin"] / (4 * pi * r**2)
                                netflux += flux

                        abc[i]["Obj_Temp"] = EquiTemp(abc[i]["Obj_Comp"], netflux)
