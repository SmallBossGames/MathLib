package smallBossMathLib.differentialEquations

interface IDifferentialEquation {
    fun getDifferential(t: Double, inY: DoubleArray, outY: DoubleArray)
}