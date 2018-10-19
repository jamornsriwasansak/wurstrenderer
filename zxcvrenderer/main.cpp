//#include "ext/tigr.h"

#include <iostream>
#include <conio.h>
#include "common/zxcv.h"

/*int main(int argc, char *argv[])
{
	Tigr *screen = tigrWindow(320, 240, "Hello", 0);
	while (!tigrClosed(screen))
	{
		tigrClear(screen, tigrRGB(0x80, 0x90, 0xa0));
		tigrPrint(screen, tfont, 120, 110, tigrRGB(0xff, 0xff, 0xff), "Test test test test test.");
		tigrUpdate(screen);
	}
	tigrFree(screen);
	return 0;
}*/


int main()
{
	Fimage result(1024, 512);
	
	result.forEachPixel([](const Vec2 & uv, Vec3 * color) {
		*color = uv[0];
	});

	Fimage::Save(result, "test.pfm");
	_getch();

	return 0;
}
