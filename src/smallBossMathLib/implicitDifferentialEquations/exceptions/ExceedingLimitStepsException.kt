package smallBossMathLib.implicitDifferentialEquations.exceptions

class ExceedingLimitStepsException : Exception() {
    override val message: String?
        get() = "Maximum steps exceeded"
}