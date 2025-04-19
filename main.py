"""

"""


import numpy as np
import math


def dTdx(coeff: float, 
         matrixSize: int, 
         matrix2dSize : int, 
         plateDisc: dict, 
         zConfines: tuple) -> np.ndarray:
    """
    This function calculates the dT/dx coefficient matrix, for the heat conduction equation.

    Args:
        coeff (float):      The coefficient for heat transfer that is thermal diffusivity over dx.
        matrixSize (int):   Size of the xyz plane, i.e. the product of the discretizations in the x, y and z direction.
        matrix2dSize (int): Size of the xy plane, i.e. the product of the discretizations in the x and y direction.
        plateDisc (dict):   A dictionary that holds how the plate is discretized in the different directions.
        zConfines (tuple):  Holds the indexes of where the plate is from and to i.e. (start,end)

    Returns:
        np.ndarray: The dT/dx coefficient matrix for the heat conduction equation
    """
    
    #The empty coefficient matrix
    coeffMatrix = np.zeros((matrixSize,matrixSize))
    
    #An empty coefficient matrix, that is as long, as the matrix is wide
    xCoeffVector = np.zeros(matrixSize)
    
    #The different coefficents lined up, with how the different temperatures affect each volume
    xCoeffVector[:3] = [coeff, -2*coeff, coeff]
    #Rolling it back to match with the centered difference method
    xCoeffVector = np.roll(xCoeffVector, -1)
    
    #The range in which is the plate and not the air
    for row in range(matrix2dSize * zConfines[0], matrix2dSize * zConfines[1]):
        if row % plateDisc["x"] == 0:
            #Forward difference method
            coeffMatrix[row] = np.roll(xCoeffVector, row + 1)
        elif row % plateDisc["x"] == plateDisc["x"] - 1:
            #Backward difference method
            coeffMatrix[row] = np.roll(xCoeffVector, row - 1)
        else:
            #Centered difference method
            coeffMatrix[row] = np.roll(xCoeffVector, row)
            
    return coeffMatrix
            
    
def dtdy(coeff: float, 
         matrixSize: int, 
         matrix2dSize: int, 
         plateDisc: dict, 
         zConfines: tuple) -> np.ndarray:
    """
    This function calculates the dT/dy coefficient matrix, for the heat conduction equation.

    Args:
        coeff (float):      The coefficient for heat transfer that is thermal diffusivity over dy.
        matrixSize (int):   Size of the xyz plane, i.e. the product of the discretizations in the x, y and z direction.
        matrix2dSize (int): Size of the xy plane, i.e. the product of the discretizations in the x and y direction.
        plateDisc (dict):   A dictionary that holds how the plate is discretized in the different directions.
        zConfines (tuple):  Holds the indexes of where the plate is from and to i.e. (start,end)

    Returns:
        np.ndarray: The dT/dy coefficient matrix for the heat conduction equation
    """
    
    #The empty coefficient matrix
    coeffMatrix = np.zeros((matrixSize,matrixSize))

    
    #dT/dy
    yCoeffVector = np.zeros(matrixSize)
    
    #The different coefficents lined up, with how the different temperatures affect each volume
    yCoeffVector[-plateDisc["x"]] = coeff
    yCoeffVector[0] = -2 * coeff
    yCoeffVector[plateDisc["x"]] = coeff
    
    #The range in which is the plate and not the air
    for row in range(matrix2dSize * zConfines[0], matrix2dSize * zConfines[1]):
        yLevel = row % matrix2dSize
        if yLevel < plateDisc["x"]:
            #Forward difference method
            coeffMatrix[row] += np.roll(yCoeffVector, row + plateDisc["x"])
        elif yLevel >= matrixSize - plateDisc["x"]:
            #Backward difference method
            coeffMatrix[row] += np.roll(yCoeffVector, row - plateDisc["x"])
        else:
            #Centered difference method
            coeffMatrix[row] += np.roll(yCoeffVector, row)
    
    return coeffMatrix

    
def dtdz(coeff: float, 
         matrixSize: int,
         matrix2dSize: int, 
         zConfines: tuple) -> np.ndarray:
    """
    This function calculates the dT/dz coefficient matrix, for the heat conduction equation.

    Args:
        coeff (float):      The coefficient for heat transfer that is thermal diffusivity over dz.
        matrixSize (int):   Size of the xyz plane, i.e. the product of the discretizations in the x, y and z direction.
        matrix2dSize (int): Size of the xy plane, i.e. the product of the discretizations in the x and y direction.
        zConfines (tuple):  Holds the indexes of where the plate is from and to i.e. (start,end)

    Returns:
        np.ndarray: The dT/dz coefficient matrix for the heat conduction equation
    """
    
    #The empty coefficient matrix
    coeffMatrix = np.zeros((matrixSize,matrixSize))

    
    #dT/dz
    zCoeffVector = np.zeros(matrixSize)
    
    #The different coefficents lined up, with how the different temperatures affect each volume
    zCoeffVector[-matrix2dSize] = coeff
    zCoeffVector[0] = -2 * coeff
    zCoeffVector[matrix2dSize] = coeff
    
    #The range in which is the plate and not the air
    for row in range(matrix2dSize * zConfines[0], matrix2dSize * zConfines[1]):
        if row % matrix2dSize == 0:
            #Forward difference method
            coeffMatrix[row] = np.roll(zCoeffVector, row + matrix2dSize)
        elif row % matrix2dSize == matrix2dSize - 1:
            #Backward difference method
            coeffMatrix[row] = np.roll(zCoeffVector, row - matrix2dSize)
        else:
            #Centered difference method
            coeffMatrix[row] = np.roll(zCoeffVector, row)

    return coeffMatrix


def heatConvCoeff(h: float, 
                  plateConductivity: float,
                  plateDensity: float, 
                  airDensity: float,
                  airConductivity: float, 
                  area: float, 
                  matrixSize: int,
                  matrix2dSize: int,
                  dx: float, dy: float, dz: float, airdz: float) -> np.ndarray:
    """
    Calculates the heat convection coefficient matrix, for the outer matrices,
    of a three dimensional heat matrix.

    Args:
        h (float):          The Convection Heat Transfer Coefficient
        k (float):          The thermal conductivity of the plate
        density (float):    Density of the plate
        area (float):       the area of the plate exposed to the air
        matrixSize (int):   Size of the xyz plane, i.e. the product of the discretizations in the x, y and z direction.
        matrix2dSize (int): Size of the xy plane, i.e. the product of the discretizations in the x and y direction.
        dx (float):         The length of the volumes in the x direction
        dy (float):         -||- y direction
        dz (float):         -||- x direction

    Returns:
        np.ndarray:         A heat convection coefficient matrix, specifically for dT/dt
    """
    
    #The empty coefficient matrix
    coeffMatrix = np.zeros((matrixSize,matrixSize))
    
    #Calculates the coefficient so that it equals the same as heat conduction
    #i.e. dT/dy or C° per second, since if not it will be watts or joules per second
    #And will not be able to connect them
    coeff = h * area
    
    #Calculates the heat that goes in or out of a point, if it is coordinately below the 
    #medium that it is "convectioning" with
    upwardsCoeffVector = np.zeros(matrixSize)
    upwardsCoeffVector[0] = coeff 
    upwardsCoeffVector[matrix2dSize] = -coeff
    
    #Now when it is above
    downwardCoeffVector = np.zeros(matrixSize)
    downwardCoeffVector[0] = -coeff
    downwardCoeffVector[matrix2dSize] = coeff
    
    #factor to convert from W or J/s to C°/s
    airJsCs = 1/(airDensity * dx * dy * airdz * airConductivity)
    plateJsCs = 1/(plateDensity * dx * dy * dz * plateConductivity)
    
    #Range is for the two singular matrices that hold the air temperatures
    #And calculating them both at once
    for row in range(matrix2dSize):
        #air
        coeffMatrix[row] = np.roll(upwardsCoeffVector, row) * airJsCs
        #plate
        coeffMatrix[row + matrix2dSize] = np.roll(downwardCoeffVector, row) * plateJsCs
        coeffMatrix[-(row + 1 + matrix2dSize)] = np.roll(upwardsCoeffVector, -(row + 1 + matrix2dSize)) * plateJsCs
        #air
        coeffMatrix[-(row + 1)] = np.roll(downwardCoeffVector, -(row + 1 + matrix2dSize)) * airJsCs
        
    return coeffMatrix





class CFPHE_plate:
    """
    
    """
    
    def __init__(plate,
                 startTempMatrix: np.ndarray,
                 coldAirstream: np.ndarray,
                 hotAirstream: np.ndarray,
                 airSpeed: float,           
                 simTime: float,
                 plateThermalDiffusivity: float,
                 plateThermalConductivity: float,
                 plateDensity: float,
                 convectionHeatTransferCoeff: float,
                 plateDimensions:dict,
                 ):
        """
        Args:
            startTempMatrix (np.ndarray):        The starting temperature distribution of the plate and air.
                                                 Where air is startTempMatrix[0] and startTempMatrix[-1]
            coldAirstream (np.ndarray):          An array that holds the temperature of the incoming cold air, each seperated by dt.
            hotAirstream (np.ndarray):           -||- hot air -||-
            airSpeed (float):                    The speed of the incoming air. (m/s)
            simTime (float):                     How many seconds the simulation should go over.
            plateThermalDiffusivity (float):     The thermal diffusivity of the plate.
            plateThermalConductivity (float):    The thermal conductivity of the plate.
            plateDensity (float):                The density of the plate
            convectionHeatTransferCoeff (float): The Convection heat transfer coefficient between the plate and air
            plateDimensions (dict):              A dict holding the dimensions of the plate in the cardinal directions
        """
        
        #The discretizations of the plate
        plate.disc = {"x": len(startTempMatrix[0,0]),
                      "y": len(startTempMatrix[0]),
                      "z": len(startTempMatrix)}
        
        plate.dim = plateDimensions
        
        plate.diffs = {"dx": plate.dim["x"] / plate.disc["x"],
                       "dy": plate.dim["y"] / plate.disc["y"],
                       "dz": plate.dim["z"] / plate.disc["z"]}
        #The amount of time it takes to get from one volume to another
        plate.diffs["dt"] = min(airSpeed / plate.diffs["dx"], airSpeed / plate.diffs["dy"])
        
        #How many dt's we need to run over the specified time
        #math.ceil() always rounds up to the nearest integer
        timeIntervals = math.ceil(simTime / plate.diffs["dt"])
        
        #Makes the matrix that holds all heat distributions at any time
        plate.heatMatrix = np.array((timeIntervals,
                                     plate.disc["z"],
                                     plate.disc["y"],
                                     plate.disc["x"]))
        #Initializing it
        plate.heatMatrix[0] = startTempMatrix
        
        #Plate properties
        plate.alpha =   plateThermalDiffusivity
        plate.k =       plateThermalConductivity
        plate.h =       convectionHeatTransferCoeff
        plate.density = plateDensity
        
        #Air properties
        plate.coldAir = coldAirstream
        plate.hotAir =  hotAirstream
        
        
        
        
        
        
        
    