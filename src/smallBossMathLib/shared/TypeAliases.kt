package smallBossMathLib.shared

typealias StationaryODE = (y: DoubleArray, outF: DoubleArray) -> Unit
typealias NonStationaryODE = (t:Double, y: DoubleArray, outF: DoubleArray) -> Unit

