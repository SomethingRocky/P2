"""

"""

import numpy as np

def xCondCoeff(coeff: float, matrixSize: int, plateDisc: dict, zConfines: tuple) -> np.ndarray:
    
    coeffMatrix = np.array((plateDisc["z"], plateDisc["y"], plateDisc["x"]))
    
    xCoeffVector = np.array((plateDisc["x"]))
    
    xCoeffVector[:3] = [coeff, -2*coeff, coeff]
    #Rolling it back to match with the centered difference method
    xCoeffVector = np.roll(xCoeffVector, -1)
    for yIndex, yVec in enumerate(coeffMatrix[zConfines[0]: zConfines[1] + 1]): 
        pass
    
def yCondCoeff(coeff: float, matrixSize: int, plateDisc: dict, zConfines: tuple) -> np.ndarray:
    
    coeffMatrix = np.array((plateDisc["z"], plateDisc["y"], plateDisc["x"]))
    
    
    
def zCondCoeff(coeff: float, matrixSize: int, plateDisc: dict, zConfines: tuple) -> np.ndarray:
    
    coeffMatrix = np.array((plateDisc["z"], plateDisc["y"], plateDisc["x"]))



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