from Numeric import *
import LinearAlgebra as linalg

_typeMap = {
    float: Float,
}

NumericArray = array
def array(data, dataType):
    return NumericArray(data, _typeMap[dataType])
