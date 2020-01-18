#pragma once

const int width = 600;
const int height = 600;

void draw_pixel(std::vector<uint8_t>&
	v, int _width, int x, int y, uint8_t r, uint8_t g, uint8_t b, uint8_t a)
{
	if (0 <= x && x < width && 0 <= y && y < height)
	{
		v[y * _width * 4 + x * 4] = r;
		v[y * _width * 4 + x * 4 + 1] = g;
		v[y * _width * 4 + x * 4 + 2] = b;
		v[y * _width * 4 + x * 4 + 3] = a;
	}
}

void draw_point(std::vector<uint8_t>& v, int _width, int ra, int x, int y,
	uint8_t r, uint8_t g, uint8_t b, uint8_t a)
{
	for (int i = -ra; i < ra; i++)
	{
		for (int j = -ra; j < ra; j++)
		{
			if (i*i + j * j <= ra * ra)
				draw_pixel(v, _width, x + i, y + j, r, g, b, a);
		}
	}
}

// 交换整数 a 、b 的值
inline void swap_int(int *a, int *b) {
	*a ^= *b;
	*b ^= *a;
	*a ^= *b;
}

// Bresenham's line algorithm
void draw_line(std::vector<uint8_t>& v, int _width, int x1, int y1, int x2, int y2, uint8_t r, uint8_t g, uint8_t b, uint8_t a)
{
	// 参数 c 为颜色值
	int dx = abs(x2 - x1),
		dy = abs(y2 - y1),
		yy = 0;

	if (dx < dy) {
		yy = 1;
		swap_int(&x1, &y1);
		swap_int(&x2, &y2);
		swap_int(&dx, &dy);
	}

	int ix = (x2 - x1) > 0 ? 1 : -1,
		iy = (y2 - y1) > 0 ? 1 : -1,
		cx = x1,
		cy = y1,
		n2dy = dy * 2,
		n2dydx = (dy - dx) * 2,
		d = dy * 2 - dx;

	if (yy) { // 如果直线与 x 轴的夹角大于 45 度
		while (cx != x2) {
			if (d < 0) {
				d += n2dy;
			}
			else {
				cy += iy;
				d += n2dydx;
			}
			draw_point(v, _width, 4, cy, cx, r, g, b, a);
			cx += ix;
		}
	}
	else { // 如果直线与 x 轴的夹角小于 45 度
		while (cx != x2) {
			if (d < 0) {
				d += n2dy;
			}
			else {
				cy += iy;
				d += n2dydx;
			}
			draw_point(v, _width, 4, cx, cy, r, g, b, a);
			cx += ix;
		}
	}
}

