import numpy as np
import matplotlib

import matplotlib.pyplot as plt

class eField:
    # Initialize
    def __init__(self, spacing: np.float32):
        self.spacing = spacing
        self.x = np.arange(0,20,spacing)
        self.field = np.empty_like(self.x)
        self.A = 10
        self.a = 2
        self.k = 0.4
        self.diffFlag = False

    # Function for applying main electric field equation
    def applyField(self, xArr: np.ndarray | None = []):
        if len(xArr) == 0:
            self.field = self.A * np.cos((np.pi*self.x)/self.a)*np.exp(-self.k*self.x)
        else:
            return self.A * np.cos((np.pi*xArr)/self.a)*np.exp(-self.k*xArr)

    # Plot the current electric field
    def plotField(self):
        plt.plot(self.x, self.field)
        plt.xlabel("Position (m)")
        plt.ylabel("Electric Field (V/m)")
        plt.title("Electric Field (V/m) vs Position (m)")
        plt.show()

    # Take derivatives of the current field
    def ddx(self):
            # self.forwardDiffField = self.field[:-1]-self.field[1:]
            self.forwardDiffField = np.diff(self.field) * 1 / self.spacing
            # self.centerDiffField = self.field[:-2]-self.field[2:]
            self.centerDiffField = np.gradient(self.field) * 1 / self.spacing
            self.analyticddxField = -(np.pi*self.A*np.exp(-self.k*self.x)*np.sin((np.pi*self.x)/self.a))/self.a-self.A*self.k*np.exp(-self.k*self.x)*np.cos((np.pi*self.x)/self.a)

    # Plot all 3 derivatives
    def plotDerivatives(self):
        plt.plot(self.x[:-1],self.forwardDiffField, c = 'r', label = "Forward Difference")
        plt.plot(self.x,self.centerDiffField, c = 'b', label = "Central Difference")
        plt.plot(self.x,self.analyticddxField, c = 'g', label = "Analytical Derivative")
        plt.xlabel("Position (m)")
        plt.ylabel("Derivative of Electric Field (V/m)")
        plt.title("Derivative of Electric Field (V/m) vs Position (m)")
        plt.legend()
        plt.show()

    # Plot errors
    def plotErrors(self):
        plt.semilogy(self.x[:-1],np.abs(self.forwardDiffField-self.analyticddxField[:-1]), c = 'r', label = "Forward Difference Error")
        plt.semilogy(self.x[1:-1],np.abs(self.centerDiffField-self.analyticddxField)[1:-1], c = 'b', label = "Central Difference Error")
        plt.xlabel("Position (m)")
        plt.ylabel("Difference (V/m)")
        plt.title("Difference (V/m) vs Position (m)")
        plt.legend()
        plt.show()

    # Run the calculations for part F
    def partF(self):
        log1 = []
        log2 = []
        for i in range(15):
            spacing = 10**(-(i+1))
            x = np.array([1-spacing,1,1+spacing])
            print(x)
            print(len(x))
            field = self.applyField(x)
            diffField = np.diff(field) * 1 / spacing
            ddxField = -(np.pi*self.A*np.exp(-self.k*x)*np.sin((np.pi*x)/self.a))/self.a-self.A*self.k*np.exp(-self.k*x)*np.cos((np.pi*x)/self.a)
            log2.append(np.abs(diffField[1] - ddxField[1]))
            log1.append(spacing)
            print(f"on loop {i}")
        log1 = np.flip(log1)
        log2 = np.flip(log2)
        plt.loglog(log1,log2, c = 'g')
        plt.xlabel("Delta x")
        plt.ylabel("Error (V/m)")
        plt.title("Error (V/m) vs Delta x")
        plt.show()


    # Carry out all functions relevant to part one of assignment
    def partOne(self):
        self.applyField()
        self.plotField()
        return


    # Carry out all functions relevant to part two of assignment
    def partTwo(self):
        self.ddx()
        self.plotDerivatives()
        self.plotErrors() # The central difference method is more accurate
        return


def part3():
    dat = np.loadtxt("/home/evilla/Downloads/fismflux20020824.dat")
    x = dat[:,4]*60+dat[:,3]*3600 # x axis in real time
    plt.plot(x,dat[:,55])
    plt.xlabel("Time (sec)")
    plt.ylabel("EUV Spectrum (W/m^2/nm)")
    plt.title("EUV Spectrum (W/m^2/nm) vs Time (sec)")
    plt.tight_layout()
    plt.show()






if __name__ == "__main__":
    field1 = eField(0.1)
    field2 = eField(0.05)
    field1.partOne()
    field1.partTwo()
    field2.partOne()
    field2.partTwo()
    field1.partF()
    part3()
