/**
 * lif outputs in 2050 and 2080
 */

function onYearEnd()
{
	console.log("GlobalEvent: on year end: " + Globals.year);

	if (Globals.year==30 || Globals.year==80)
        // output at 10m resolution, 2m height
		Globals.gridToFile('lifc','output/lif_' + Globals.year + '.asc', 2);
}


