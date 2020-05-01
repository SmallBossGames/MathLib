package smallBossMathLib.explicitDifferentialEquations

data class ExplicitMethodStatistic(
    override var stepsCount: Int,
    override var evaluationsCount: Int,
    override var returnsCount: Int
) : IExplicitMethodStatistic