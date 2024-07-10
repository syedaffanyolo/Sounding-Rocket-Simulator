import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np

LARGEFONT = ("Verdana", 15)
# Given parameters
g = 9.81  # m/s^2, acceleration due to gravity

cd = 0.64  # drag coefficient for truncated tangent ogive nose cone
cd2 = 1.0  # drag coefficient for boost phase
cdp = 1.75  # drag coefficient for recovery under parachute phase

class Store:
    def __init__(self):
        self.param = ""  # Initialize the variable
        self.isp = 122.96
        self.Rmo = 12.137
        self.thrust = 2006
        self.Pmo = 3.34956
        self.p = 100575.7
        self.t = 26
        self.b = 0.096
        self.l = 2.3
        self.h = 0.3
        self.dp = 1.5
        self.smax = 0

    def set_param(self, new_string, isp, Rmo, thrust, Pmo, p, t, b, l, h, dp):
        self.param = new_string  # Set the variable to a new string
        self.isp = isp
        self.Rmo = Rmo
        self.thrust = thrust
        self.Pmo = Pmo
        self.p = p
        self.t = t
        self.b = b
        self.l = l
        self.h = h
        self.dp = dp

    def set_smax(self, smax):
        self.smax = smax

    def get_smax(self):
        return self.smax

    def get_string(self):
        return self.param  # Access the variable
    
    def get_param(self, parameter):
        return getattr(self, parameter, None)

store = Store()

class tkinterApp(tk.Tk):
    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)
        
        self.resize_window(400, 650)
        
        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)
        
        self.frames = {}
        
        for F in (StartPage, Page1):
            frame = F(container, self)
            self.frames[F] = frame
            frame.grid(row=0, column=0, sticky="nsew")

        self.show_frame(StartPage)

    def show_frame(self, cont):
        frame = self.frames[cont]
        frame.tkraise()

    def resize_window(self, width, height):
        self.geometry(f"{width}x{height}")


class StartPage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        
        def action():
            controller.resize_window(1000, 1000)
            controller.show_frame(Page1)
            isp = self.parameters["Specific Impulse (s)"].get()
            Rmo = self.parameters["Rocket Mass (kg)"].get()
            thrust = self.parameters["Thrust (N)"].get()
            Pmo = self.parameters["Propellant Mass (kg)"].get()
            p = self.parameters["Atmosphere Pressure (inHg)"].get()
            t = self.parameters["Temperature (C)"].get()
            b = self.parameters["Airframe Diameter (m)"].get()
            l = self.parameters["Airframe Length (m)"].get()
            h = self.parameters["Nosecone Length (m)"].get()
            dp = self.parameters["Parachute Diameter (m)"].get()
            store.set_param(combo.get(), isp, Rmo, thrust, Pmo, p, t, b, l, h, dp)
            controller.frames[Page1].display_selection()
            
        self.labelHeading = ttk.Label(self, text="Enter Parameters", font=LARGEFONT)
        self.label = ttk.Label(self, text="")

        self.labelHeading.place(x=140, y=50)
        self.label.place(x=55, y=100)

        combo = ttk.Combobox(self, state="readonly", values=["Altitude vs Time", "Velocity vs Time", "Drag Force vs Time"])
        combo.grid(row=1, column=0)
        combo.place(x=90, y=500)

        self.parameters = {
            "Thrust (N)": tk.DoubleVar(),
            "Rocket Mass (kg)": tk.DoubleVar(),
            "Propellant Mass (kg)": tk.DoubleVar(),
            "Specific Impulse (s)": tk.DoubleVar(),
            "Atmosphere Pressure (inHg)": tk.DoubleVar(),
            "Temperature (C)": tk.DoubleVar(),
            "Airframe Length (m)": tk.DoubleVar(),
            "Nosecone Length (m)": tk.DoubleVar(),
            "Airframe Diameter (m)": tk.DoubleVar(),
            "Parachute Diameter (m)": tk.DoubleVar(),
        }

        row = 0
        for param, var in self.parameters.items():
            ttk.Label(self.label, text=param).grid(row=row, column=0, padx=5, pady=5, sticky=tk.W)
            ttk.Entry(self.label, textvariable=var, width=10).grid(row=row, column=1, padx=5, pady=5)
            row += 1

        button1 = ttk.Button(self, text="Submit", command=action)
        button1.place(x=140, y=550)


class Page1(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)

        self.controller = controller

        self.fig = Figure(figsize=(8, 6), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        button = ttk.Button(self, text="Back", command=self.go_back)
        button.pack(side=tk.BOTTOM)

    def go_back(self):
        self.controller.resize_window(400, 650)
        self.controller.show_frame(StartPage)

    def plot_alt(self, x, y):
        self.ax.clear()
        self.ax.plot(x, y)
        self.ax.set_xlabel('Time (s)')
        self.ax.set_ylabel('Altitude (m)')
        self.ax.set_title('Plot of Altitude vs Time')
        self.ax.grid(True)
        self.canvas.draw()

    def plot_vel(self, x, y):
        self.ax.clear()
        self.ax.plot(x, y)
        self.ax.set_xlabel('Time (s)')
        self.ax.set_ylabel('Velocity (m/s)')
        self.ax.set_title('Plot of Velocity vs Time')
        self.ax.grid(True)
        self.canvas.draw()

    def plot_drg(self, x, y):
        self.ax.clear()
        self.ax.plot(x, y)
        self.ax.set_xlabel('Time (s)')
        self.ax.set_ylabel('Drag Force (N)')
        self.ax.set_title('Plot of Drag Force vs Time')
        self.ax.grid(True)
        self.canvas.draw()

    def display_selection(self):
        selection = store.get_string()  # Get the selection from the store
        print(selection)  # Print the selection stored in the store

        # Calculated parameters
        rho = store.get_param("p") /(287.058*(store.get_param("t") + 273.15))
        c = store.get_param("isp") * g  # m/s, effective exhaust velocity
        A = np.pi * ((store.get_param("b") / 2) ** 2)  # m^2, cross-sectional area of airframe from seen from top
        A2 = (store.get_param("b") * store.get_param("l")) + (0.5 * store.get_param("b") * store.get_param("h"))
        Ap = np.pi * ((store.get_param("dp") / 2) ** 2)  # cross-section area of parachute
        mbo = store.get_param("Rmo") - store.get_param("Pmo")  # kg, burnout mass
        mdot = store.get_param("thrust") / (store.get_param("isp") * g)  # kg/s, mass burn rate
        tb = store.get_param("Pmo") / mdot  # s, time at burnout

        k = cd * A * 0.5 * rho
        k2 = cd2 * A2 * 0.5 * rho
        kp = cdp * Ap * 0.5 * rho

        # Burn phase calculations
        psi = k / mbo
        A1, A2, A3 = [], [], []

        for dt in np.arange(0, tb, 0.03):
            W = store.get_param("Rmo") - (mdot * dt)
            v = (c * np.log(store.get_param("Rmo") / W)) - (g * dt)
            D = k * v**2

            cosh_arg = ((dt / W) * np.sqrt(k * (store.get_param("thrust") - (W * g))))
            sb = (W / k) * np.log(np.cosh(cosh_arg))

            A1.append(sb)
            A2.append(v)
            A3.append(D)

        tc = np.sqrt(mbo / (g * k)) * np.arctan(v * np.sqrt(k / (mbo * g)))

        for dt1 in np.arange(0, tc, 0.03):
            vc = np.sqrt(g / psi) * np.tan(np.arctan(np.sqrt(psi / g) * v) - np.sqrt(psi * g) * dt1)
            sc = (mbo / (2 * k)) * np.log(((mbo * g) + (k * (v - vc)**2)) / (mbo * g))
            Dc = k * vc**2
            store.set_smax(sc + sb) 
            A1.append(store.get_smax())
            A2.append(vc)
            A3.append(Dc)

        # Descent phase calculations
        tt = store.get_smax() * np.sqrt(kp / (mbo * g))

        for dt2 in np.arange(0, tt, 0.03):
            exponent = np.sqrt(2 * cdp * Ap * rho * g / mbo) * dt2
            numerator = np.sqrt(2 * g * mbo) * (np.exp(exponent) - 1)
            denominator = np.sqrt(cdp * Ap * rho) * (np.exp(exponent) + 1)
            vf = numerator / denominator
            vt = np.sqrt(mbo * g / kp)
            Df = kp * vf**2
            sf = store.get_smax() - vt * dt2 + ((vt * mbo) / kp) * (np.exp(-kp * dt2 / mbo) - 1)

            if sf > 0:
                A1.append(sf)
                A2.append(vf)
                A3.append(Df)

        t = np.arange(0, len(A1) * 0.03, 0.03)
        t2 = np.arange(0, len(A2) * 0.03, 0.03)
        t3 = np.arange(0, len(A3) * 0.03, 0.03)

        if selection == "Altitude vs Time":
            self.plot_alt(t, A1)
        elif selection == "Velocity vs Time":
            self.plot_vel(t2, A2)
        elif selection == "Drag Force vs Time":
            self.plot_drg(t3, A3)


app = tkinterApp()
app.mainloop()
