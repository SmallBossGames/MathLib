package smallBossMathLib.shared

interface IImplicitMethodStepInfo {
    val time:Double
    val yValue:DoubleArray
    val isLowLimitReached:Boolean
    val isHighLimitReached: Boolean
    val stepsCount: Int
    val evaluationsCount: Int
    val jacobiEvaluationsCount: Int
    val returnsCount: Int
}