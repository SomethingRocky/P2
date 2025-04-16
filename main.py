"""

"""

import numpy as np

def dTdx(coeff: float, 
         matrixSize: int, 
         matrix2dSize : int, 
         plateDisc: dict, 
         zConfines: tuple) -> np.ndarray:
    
    
    coeffMatrix = np.zeros((matrixSize,matrixSize))
    
    #dT/dx
    xCoeffVector = np.zeros(matrixSize)
    
    xCoeffVector[:3] = [coeff, -2*coeff, coeff]
    #Rolling it back to match with the centered difference method
    xCoeffVector = np.roll(xCoeffVector, -1)
    
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
    
    coeffMatrix = np.zeros((matrixSize,matrixSize))

    
    #dT/dy
    yCoeffVector = np.zeros(matrixSize)
    

    yCoeffVector[-plateDisc["x"]] = coeff
    yCoeffVector[0] = -2 * coeff
    yCoeffVector[plateDisc["x"]] = coeff
    
 
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
         plateDisc: dict, 
         zConfines: tuple) -> np.ndarray:
    
    coeffMatrix = np.zeros((matrixSize,matrixSize))

    
    #dT/dz
    zCoeffVector = np.zeros(matrixSize)
    
    zCoeffVector[-matrix2dSize] = coeff
    zCoeffVector[0] = -2 * coeff
    zCoeffVector[matrix2dSize] = coeff
    
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
                  k: float,
                  density: float,  
                  area: float, 
                  matrixSize: int,
                  matrix2dSize: int,
                  dx: float, dy: float, dz: float):
    """"""
    coeffMatrix = np.zeros((matrixSize,matrixSize))
    
    coeff = (h * area)/(density * dx * dy * dz * k)
    
    aboveCoeffVector = np.zeros(matrixSize)
    aboveCoeffVector[0] = coeff 
    aboveCoeffVector[matrix2dSize] = -coeff
    
    belowCoeffVector = np.zeros(matrixSize)
    belowCoeffVector[0] = -coeff
    belowCoeffVector[matrix2dSize] = coeff
    
    for row in range(matrix2dSize):
        coeffMatrix[row] = np.roll(aboveCoeffVector, row)
        coeffMatrix[row + matrix2dSize] = np.roll(belowCoeffVector, row)
        coeffMatrix[-(row + 1 + matrix2dSize)] = np.roll(aboveCoeffVector, -(row + 1 + matrix2dSize))
        coeffMatrix[-(row + 1)] = np.roll(belowCoeffVector, -(row + 1 + matrix2dSize))
        
    return coeffMatrix


class CFPHE_plate:
    """
    
    """
    
    def __init__(plate,
                 startTempMatrix: np.ndarray,
                 simTime: float,
                 ):
        pass
    