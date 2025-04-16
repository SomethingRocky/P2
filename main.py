"""

"""

import numpy as np

def xCondCoeff(coeff: float, matrixSize: int, plateDisc: dict, zConfines: tuple) -> np.ndarray:
    
    coeffMatrix = np.zeros((matrixSize,matrixSize))
    Matrix2dSize = plateDisc["x"] * plateDisc["y"]
    
    #dT/dx
    xCoeffVector = np.zeros(matrixSize)
    
    xCoeffVector[:3] = [coeff, -2*coeff, coeff]
    #Rolling it back to match with the centered difference method
    xCoeffVector = np.roll(xCoeffVector, -1)
    
    for row in range(Matrix2dSize * zConfines[0], Matrix2dSize * zConfines[1]):
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
            
    
def yCondCoeff(coeff: float, matrixSize: int, plateDisc: dict, zConfines: tuple) -> np.ndarray:
    
    coeffMatrix = np.zeros((matrixSize,matrixSize))
    Matrix2dSize = plateDisc["x"] * plateDisc["y"]
    
    #dT/dy
    yCoeffVector = np.zeros(matrixSize)
    

    yCoeffVector[-plateDisc["x"]] = coeff
    yCoeffVector[0] = -2 * coeff
    yCoeffVector[plateDisc["x"]] = coeff
    
 
    for row in range(Matrix2dSize * zConfines[0], Matrix2dSize * zConfines[1]):
        yLevel = row % Matrix2dSize
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
    
    
    
def zCondCoeff(coeff: float, matrixSize: int, plateDisc: dict, zConfines: tuple) -> np.ndarray:
    
    coeffMatrix = np.zeros((matrixSize,matrixSize))
    matrix2dSize = plateDisc["x"] * plateDisc["y"]
    
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



def heatCondCoeff(α: float, plateDisc: dict, dx: float, dy: float, dz: float):
    
    matrixSize = plateDisc["x"] * plateDisc["y"] * plateDisc["z"]
    

def heatConvCoeff(h: float, area: float, plateDisc: dict, dx: float, dy: float, dz: float):
    pass


class CFPHE_plate:
    """
    
    """
    
    def __init__(plate,
                 startTempMatrix: np.ndarray,
                 simTime: float,
                 ):
        pass