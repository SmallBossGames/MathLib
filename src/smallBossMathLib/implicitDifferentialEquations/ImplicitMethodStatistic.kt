package smallBossMathLib.implicitDifferentialEquations

data class ImplicitMethodStatistic(
    override var stepsCount: Int,
    override var evaluationsCount: Int,
    override var jacobiEvaluationsCount: Int,
    override var returnsCount: Int
) : IImplicitMethodStatistic