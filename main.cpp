#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <ctime>
#include <Windows.h>
#include <MMSystem.h>

class Wavelet {
private:
	inline int I(const int y, const int x) {
		return y * size_ + x;
	}

	// level = 0; 全体
	// level = 1; (1/2)x(1/2)
	// level = 2; (1/4)x(1/4) ...
	bool wavelet2D(const double *input, double *output, const int now_size) {
		const int now_size_half = now_size / 2;

		if (now_size <= 1 || size_ < now_size)
			return false;

		std::vector<double> tmp(now_size);

		// 横方向
		for (int y = 0; y < now_size; ++y) {
			for (int x = 0; x < now_size_half; ++x) {
				const double v0 = input[I(y, 2 * x)];
				const double v1 = input[I(y, 2 * x + 1)];
				output[I(y, x)] = (v0 + v1) / 2;
				output[I(y, x + now_size_half)] = (v0 - v1) / 2;
			}
		}
	
		// 縦方向
		for (int x = 0; x < now_size; ++x) {
			for (int y = 0; y < now_size; ++y)
				tmp[y] = output[I(y, x)];
			for (int y = 0; y < now_size_half; ++y) {
				const double v0 = tmp[2 * y];
				const double v1 = tmp[2 * y + 1];
				output[I(y, x)] = (v0 + v1) / 2;
				output[I(y + now_size_half, x)] = (v0 - v1) / 2;
			}
		}
		return true;
	}

	bool wavelet2Dreverse(const double *input, double *output, const int now_size) {
		const int now_size_half = now_size / 2;
		
		if (now_size <= 1 || size_ < now_size)
			return false;

		std::vector<double> tmp(now_size);
	
	
		// 縦方向
		for (int x = 0; x < now_size; ++x) {
			for (int y = 0; y < now_size_half; ++y) {
				const double v0 = input[I(y, x)];
				const double v1 = input[I(y + now_size_half, x)];
				output[I(2 * y, x)] = v0 + v1;
				output[I(2 * y + 1, x)] = v0 - v1;
			}
		}

		// 横方向
		for (int y = 0; y < now_size; ++y) {
			for (int x = 0; x < now_size; ++x)
				tmp[x] = output[I(y, x)];
			for (int x = 0; x < now_size_half; ++x) {
				const double v0 = tmp[x];
				const double v1 = tmp[x + now_size_half];
				output[I(y, 2 * x)] = v0 + v1;
				output[I(y, 2 * x + 1)] = v0 - v1;
			}
		}
		return true;
	}

	int size_;
	std::vector<double> image_;

	std::vector<double> transformed_;
	std::vector<double> reversed_;


	struct Coefficient {
		unsigned char index_;
		double val_;
		Coefficient(const unsigned char index, double val) :
		index_(index), val_(val) {}
	};

	std::vector<Coefficient> coefficients_; // 非ゼロWavelet係数を一列に並べたもの
	std::vector<int> num_coefficient_; // ブロックごとの非ゼロWavelet係数の数
	int num_block_; // ブロック数
public:
	Wavelet(const std::string &filename) {
		size_ = 0;
		num_block_ = 0;

		FILE *fp = fopen(filename.c_str(), "rb");
		if (fp != NULL) {
			int ch;
			std::vector<unsigned char> rawdata;
			while ((ch = fgetc(fp)) != EOF) {
				rawdata.push_back(ch);
			}
			fclose(fp);
			image_.resize(rawdata.size() / 3);
			for (int i = 0; i < rawdata.size() / 3; ++i)
				image_[i] = (double)rawdata[i * 3] /* / 255.0 */;
			size_ = (int)sqrt((double)image_.size());
		}
	}

	void save(const std::string &filename, const std::vector<double> data) {
		std::vector<unsigned char> final_image(3 * data.size());
		for (int i = 0; i < data.size(); i ++) {
			unsigned char uc = 0;
			const double d = data[i] * 255.0;
			if (d < 0.0)
				uc = 0;
			else if (d >= 255.0)
				uc = 255;
			else
				uc = (unsigned char)d;

			final_image[i*3] = uc;
			final_image[i*3+1] = uc;
			final_image[i*3+2] = uc;
		}
		
		FILE *fp = fopen(filename.c_str(), "wb");
		if (fp != NULL) {
			fwrite(&final_image[0], sizeof(unsigned char), final_image.size(), fp);
			fclose(fp);
		}
	}

	void save_reversed(const std::string &filename) {
		save(filename, reversed_);
	}

	void cut(const double alpha) {
		std::vector<double> workspace(transformed_.size());
		for (int i = 0; i < workspace.size(); i ++)
			workspace[i] = abs(transformed_[i]);
		std::sort(workspace.begin(), workspace.end());
		const double threashold = workspace[(unsigned int)(alpha * workspace.size())];

		for (int i = 0; i < transformed_.size(); i ++) {
			if (abs(transformed_[i]) <= threashold) {
				transformed_[i] = 0.0;
			}
		}
	}

	void transform() {
		std::vector<double> output = image_;
		std::vector<double> output_tmp = image_;

		int now_size = size_;
		for (int i = 0; ; ++i) {
			if (!wavelet2D(&output[0], &output_tmp[0], now_size))
				break;
			now_size /= 2;
			output = output_tmp;
		}

		transformed_ = output;
	}


	void unnormalize() {
		int size = 1;
		for (int level = 0; ; ++level) {
			for (int i = 0; i < size; ++i) {
				for (int j = 0; j < size; ++j) {
					int x10, y10;
					int x01, y01;
					int x11, y11;

					raw_index(level, i, j, 0x10, &x10, &y10);
					raw_index(level, i, j, 0x01, &x01, &y01);
					raw_index(level, i, j, 0x11, &x11, &y11);

					transformed_[I(y01, x01)] *= (1 << level);
					transformed_[I(y10, x10)] *= (1 << level);
					transformed_[I(y11, x11)] *= (1 << level);
				}
			}
			size *= 2;
			if (size == size_)
				break;
		}
	}
	void normalize() {
		int size = 1;
		for (int level = 0; ; ++level) {
			for (int i = 0; i < size; ++i) {
				for (int j = 0; j < size; ++j) {
					int x10, y10;
					int x01, y01;
					int x11, y11;

					raw_index(level, i, j, 0x10, &x10, &y10);
					raw_index(level, i, j, 0x01, &x01, &y01);
					raw_index(level, i, j, 0x11, &x11, &y11);

					const double w = pow(0.5, level);
					transformed_[I(y01, x01)] *= w;
					transformed_[I(y10, x10)] *= w;
					transformed_[I(y11, x11)] *= w;
				}
			}
			size *= 2;
			if (size == size_)
				break;
		}
	}

	void sparsing() {
		num_block_ = transformed_.size() / 256;
		int now_block = 0;
		int counter = 0;
		for (int i = 0; i < transformed_.size(); i ++) {
			if (fabs(transformed_[i]) > 0.0) { // 非ゼロ
				coefficients_.push_back(Coefficient(i % 256, transformed_[i]));
				counter ++;
			}

			if (i % 256 == 255) {
				num_coefficient_.push_back(counter);
				counter = 0;
			}
		}
		/*
		std::cout << num_block_ << std::endl;
		std::cout << num_coefficient_.size() << std::endl;
		std::cout << coefficients_.size() << std::endl;
		*/
	}

	// M
	// 00 10
	// 01 11


	// 10
	// □■
	// □■

	// 01
	// □□
	// ■■
	
	// 11
	// □■
	// ■□
	void raw_index(const int level, const int offset_i, const int offset_j, const int M, int *x, int *y) {
		if (M == 0x00) {
			*x = *y = 0;
			return;
		}

		const int offset = 1 << level;
		if (M == 0x10) {
			*x = offset + offset_i;
			*y = offset_j;
		} else if (M == 0x01) {
			*x = offset_i;
			*y = offset + offset_j;
		} else if (M == 0x11) {
			*x = offset + offset_i;
			*y = offset + offset_j;
		}
	}

	int sign_of_quadrant(const int M, const int x, const int y) {
		const int table[4][2][2] = {
			{
				{1, 1},
				{1, 1}
			},
			{
				{1, 1},
				{-1, -1}
			},
			{
				{1, -1},
				{1, -1}
			},
			{
				{1, -1},
				{-1, 1}
			},
		};

		switch (M) {
		case 0x00:
			return table[0][y][x];
		case 0x01:
			return table[1][y][x];
		case 0x10:
			return table[2][y][x];
		case 0x11:
			return table[3][y][x];
		}
	}

	std::vector<double> psum_table;
	std::vector<double> csum_table;

	double psum(Wavelet &w, const int l, const int i, const int j) {
		const int idx = j * size_ + i + (1 << l);
		if (w.psum_table[idx] < 1e255) {
			return w.psum_table[idx];
		}

		if (l == 0 && i == 0 && j == 0)
			return w.transformed_[0];

		const int ol = l - 1;
		const int oi = i / 2;
		const int oj = j / 2;

		const int qx = i - 2 * oi;
		const int qy = j - 2 * oj;
		
		int x10, y10;
		int x01, y01;
		int x11, y11;

		raw_index(ol, oi, oj, 0x10, &x10, &y10);
		raw_index(ol, oi, oj, 0x01, &x01, &y01);
		raw_index(ol, oi, oj, 0x11, &x11, &y11);


		const double ret = psum(w, ol, oi, oj) + (1 << ol) * (
			w.transformed_[y01 * size_ + x01] * sign_of_quadrant(0x01, qx, qy) +
			w.transformed_[y10 * size_ + x10] * sign_of_quadrant(0x10, qx, qy) +
			w.transformed_[y11 * size_ + x11] * sign_of_quadrant(0x11, qx, qy));

		w.psum_table[idx] = ret;
		return ret;
	}

	
	void up(const int M, std::vector<double> &result, Wavelet &w1, Wavelet &w2, const int l, const int i, const int j) {
		int ol = l - 1;
		int oi = i / 2;
		int oj = j / 2;

		int ooi = i;
		int ooj = j;
		
		int x, y;
		raw_index(l, i, j, M, &x, &y);

		const double value = w1.transformed_[y * size_ + x] * w2.transformed_[y * size_ + x];
		result[0] += value;

		for (;;) {
			if (ol < 0)
				break;

			const int qx = ooi - 2 * oi;
			const int qy = ooj - 2 * oj;
		
			int x10, y10;
			int x01, y01;
			int x11, y11;

			raw_index(ol, oi, oj, 0x10, &x10, &y10);
			raw_index(ol, oi, oj, 0x01, &x01, &y01);
			raw_index(ol, oi, oj, 0x11, &x11, &y11);

			result[y01 * size_ + x01] += (1 << ol) * value * sign_of_quadrant(0x01, qx, qy);
			result[y10 * size_ + x10] += (1 << ol) * value * sign_of_quadrant(0x10, qx, qy);
			result[y11 * size_ + x11] += (1 << ol) * value * sign_of_quadrant(0x11, qx, qy);

			ol --;
			ooi = oi;
			ooj = oj;
			oi /= 2;
			oj /= 2;
		}

	}
	
	double csum2(Wavelet &w1, Wavelet &w2, const int l, const int i, const int j) {
		const int idx = j * size_ + i + (1 << l);
		if (w1.csum_table[idx] < 1e255) {
			return w1.csum_table[idx];
		}
		const int cl = l + 1;
		double tmp = 0.0;
		if ((1 << cl) >= size_)
			return 0.0;

		const int offset_i = i * 2;
		const int offset_j = j * 2;
		for (int ci = 0; ci < 2; ci ++) {
			for (int cj = 0; cj < 2; cj ++) {
				int x10, y10;
				int x01, y01;
				int x11, y11;

				raw_index(cl, offset_i + ci, offset_j + cj, 0x10, &x10, &y10);
				raw_index(cl, offset_i + ci, offset_j + cj, 0x01, &x01, &y01);
				raw_index(cl, offset_i + ci, offset_j + cj, 0x11, &x11, &y11);
				
				tmp +=
					w1.transformed_[y01 * size_ + x01] * w2.transformed_[y01 * size_ + x01] +
					w1.transformed_[y10 * size_ + x10] * w2.transformed_[y10 * size_ + x10] +
					w1.transformed_[y11 * size_ + x11] * w2.transformed_[y11 * size_ + x11];

				tmp += csum2(w1, w2, cl, offset_i + ci, offset_j + cj);
			}
		}
		w1.csum_table[idx] = tmp;
		return tmp;
	}
	double csum(const int M, Wavelet &w1, Wavelet &w2, const int l, const int i, const int j) {
		const int cl = l + 1;
		double tmp = 0.0;
		if ((1 << cl) >= size_)
			return 0.0;

		const int offset_i = i * 2;
		const int offset_j = j * 2;
		for (int ci = 0; ci < 2; ci ++) {
			for (int cj = 0; cj < 2; cj ++) {
				int x10, y10;
				int x01, y01;
				int x11, y11;

				raw_index(cl, offset_i + ci, offset_j + cj, 0x10, &x10, &y10);
				raw_index(cl, offset_i + ci, offset_j + cj, 0x01, &x01, &y01);
				raw_index(cl, offset_i + ci, offset_j + cj, 0x11, &x11, &y11);
				
				tmp += (1 << l) * (
					w1.transformed_[y01 * size_ + x01] * w2.transformed_[y01 * size_ + x01] +
					w1.transformed_[y10 * size_ + x10] * w2.transformed_[y10 * size_ + x10] +
					w1.transformed_[y11 * size_ + x11] * w2.transformed_[y11 * size_ + x11] +
					csum2(w1, w2, cl, offset_i + ci, offset_j + cj)) * sign_of_quadrant(M, ci, cj);
			}
		}
		return tmp;
	}

	int nlz2(unsigned int x)
	{
		unsigned int y;
		int n = 32;
		y = x >> 16; if (y != 0){ n = n - 16 ; x = y; }
		y = x >>  8; if (y != 0){ n = n -  8 ; x = y; }
		y = x >>  4; if (y != 0){ n = n -  4 ; x = y; }
		y = x >>  2; if (y != 0){ n = n -  2 ; x = y; }
		y = x >>  1; if (y != 0){ return n - 2; }
		return n - x;
	}

	void wavelet_index(const int x, const int y, int *M, int *level, int *i, int *j) {
		const int nlz_x = 32 - nlz2(x);
		const int nlz_y = 32 - nlz2(y);
//				std::cout << nlz_x << " " << nlz_y << std::endl;

		if (nlz_x <= 1 && nlz_y <= 1)
			*level = 0;
		else
			*level = max(nlz_x, nlz_y) - 1;

		const int offset = 1 << *level;

		*M = 0;
		if (x >= offset) {
			*M = *M | 0x10;
			*i = x - offset;
		} else {
			*i = x;
		}
		if (y >= offset) {
			*M = *M | 0x01;
			*j = y - offset;
		} else {
			*j = y;
		}
	}
	void product_sublinear(Wavelet &w) {
		// sumテーブル作る
		//cnt = 0;
		psum_table.resize(transformed_.size());
		w.psum_table.resize(transformed_.size());
		/*
		csum_table.resize(transformed_.size());
		w.csum_table.resize(transformed_.size());
		*/

		const int tmp = psum_table.size();
		for (int i = 0; i < tmp; ++i) {
			psum_table[i] = 1e255;
			w.psum_table[i] = 1e255;
			/*
			csum_table[i] = 1e255;
			w.csum_table[i] = 1e255;
			*/
		}

		std::vector<double> result(transformed_.size());
		// case 1
		// up()内部で処理
		//result[0] += transformed_[0] * w.transformed_[0];
		// case 2
		// スパース
		int lg = 0;
		int count = 0;
		for (int ib = 0; ib < num_block_; ++ib) {
			if (num_coefficient_[ib] == 0)
				continue;
			for (int ic = 0; ic < num_coefficient_[ib]; ++ic) {
				lg ++;
				const int index = coefficients_[count + ic].index_ + 256 * ib;
				const int x = index % size_;
				const int y = index / size_;

				int M, level, i, j;

				wavelet_index(x, y, &M, &level, &i, &j);
				const int Cuvw = 1 << level;
				
				int x10, y10;
				int x01, y01;
				int x11, y11;

				raw_index(level, i, j, 0x10, &x10, &y10);
				raw_index(level, i, j, 0x01, &x01, &y01);
				raw_index(level, i, j, 0x11, &x11, &y11);

				switch (M) {
				case 0x00:
					break;

				case 0x01:
					result[y11 * size_ + x11] += Cuvw * coefficients_[count + ic].val_ * w.transformed_[y10 * size_ + x10];
					result[y10 * size_ + x10] += Cuvw * coefficients_[count + ic].val_ * w.transformed_[y11 * size_ + x11];
					break;

				case 0x10:
					result[y11 * size_ + x11] += Cuvw * coefficients_[count + ic].val_ * w.transformed_[y01 * size_ + x01];
					result[y01 * size_ + x01] += Cuvw * coefficients_[count + ic].val_ * w.transformed_[y11 * size_ + x11];
					break;

				case 0x11:
					result[y01 * size_ + x01] += Cuvw * coefficients_[count + ic].val_ * w.transformed_[y10 * size_ + x10];
					result[y10 * size_ + x10] += Cuvw * coefficients_[count + ic].val_ * w.transformed_[y01 * size_ + x01];
					break;
				}
			}
			count += num_coefficient_[ib];
		}


		/*
		// インデックスベースで
		for (int y = 0; y < size_; ++y) {
			for (int x = 0; x < size_; ++x) {
				int M, level, i, j;

				wavelet_index(x, y, &M, &level, &i, &j);
				const int Cuvw = 1 << level;
				
				int x10, y10;
				int x01, y01;
				int x11, y11;

				raw_index(level, i, j, 0x10, &x10, &y10);
				raw_index(level, i, j, 0x01, &x01, &y01);
				raw_index(level, i, j, 0x11, &x11, &y11);

				switch (M) {
				case 0x00:
					break;

				case 0x01:
					result[y11 * size_ + x11] += Cuvw * transformed_[y * size_ + x] * w.transformed_[y10 * size_ + x10];
					result[y10 * size_ + x10] += Cuvw * transformed_[y * size_ + x] * w.transformed_[y11 * size_ + x11];
					break;

				case 0x10:
					result[y11 * size_ + x11] += Cuvw * transformed_[y * size_ + x] * w.transformed_[y01 * size_ + x01];
					result[y01 * size_ + x01] += Cuvw * transformed_[y * size_ + x] * w.transformed_[y11 * size_ + x11];
					break;

				case 0x11:
					result[y01 * size_ + x01] += Cuvw * transformed_[y * size_ + x] * w.transformed_[y10 * size_ + x10];
					result[y10 * size_ + x10] += Cuvw * transformed_[y * size_ + x] * w.transformed_[y01 * size_ + x01];
					break;
				}
			}
		}*/

		// case 3
		
		// スパース
		/*
		// up()内部で処理
		result[0] += transformed_[1] * w.transformed_[1];
		result[0] += transformed_[size_] * w.transformed_[size_];
		result[0] += transformed_[size_ + 1] * w.transformed_[size_ + 1];
		result[0] += csum(0x00, w, *this, 0, 0, 0); */
		count = 0;
		for (int ib = 0; ib < num_block_; ++ib) {
			if (num_coefficient_[ib] == 0)
				continue;
			for (int ic = 0; ic < num_coefficient_[ib]; ++ic) {
				lg ++;
				const int index = coefficients_[count + ic].index_ + 256 * ib;
				const int x = index % size_;
				const int y = index / size_;

				int M, level, i, j;
				wavelet_index(x, y, &M, &level, &i, &j);
				
				int x10, y10;
				int x01, y01;
				int x11, y11;
				
				switch (M) {
				case 0x00:
					break;

				case 0x01:
					raw_index(level, i, j, 0x01, &x01, &y01);
					result[y01 * size_ + x01] += coefficients_[count + ic].val_ * psum(w, level, i, j);
					break;

				case 0x10:
					raw_index(level, i, j, 0x10, &x10, &y10);
					result[y10 * size_ + x10] += coefficients_[count + ic].val_ * psum(w, level, i, j);
					break;

				case 0x11:
					raw_index(level, i, j, 0x11, &x11, &y11);
					result[y11 * size_ + x11] += coefficients_[count + ic].val_ * psum(w, level, i, j);
					break;
				}

			}
			count += num_coefficient_[ib];
		}
		
		count = 0;
		for (int ib = 0; ib < w.num_block_; ++ib) {
			if (w.num_coefficient_[ib] == 0)
				continue;
			for (int ic = 0; ic < w.num_coefficient_[ib]; ++ic) {
				lg ++;
				const int index = w.coefficients_[count + ic].index_ + 256 * ib;
				const int x = index % size_;
				const int y = index / size_;

				int M, level, i, j;
				wavelet_index(x, y, &M, &level, &i, &j);
				
				int x10, y10;
				int x01, y01;
				int x11, y11;
				
				switch (M) {
				case 0x00:
					break;

				case 0x01:
					raw_index(level, i, j, 0x01, &x01, &y01);
					result[y01 * size_ + x01] += w.coefficients_[count + ic].val_ * psum(*this, level, i, j);
					break;

				case 0x10:
					raw_index(level, i, j, 0x10, &x10, &y10);
					result[y10 * size_ + x10] += w.coefficients_[count + ic].val_ * psum(*this, level, i, j);
					break;

				case 0x11:
					raw_index(level, i, j, 0x11, &x11, &y11);
					result[y11 * size_ + x11] += w.coefficients_[count + ic].val_ * psum(*this, level, i, j);
					break;
				}

			}
			count += w.num_coefficient_[ib];
		}
		
		count = 0;
		for (int ib = 0; ib < w.num_block_; ++ib) {
			if (w.num_coefficient_[ib] == 0)
				continue;
			for (int ic = 0; ic < w.num_coefficient_[ib]; ++ic) {
				lg ++;
				const int index = w.coefficients_[count + ic].index_ + 256 * ib;
				const int x = index % size_;
				const int y = index / size_;

				int M, level, i, j;
				wavelet_index(x, y, &M, &level, &i, &j);

				up(M, result, *this, w, level, i, j);
			}
			count += w.num_coefficient_[ib];
		}

		transformed_ = result;
	}

	void product_linear(Wavelet &w) {
		// sumテーブル作る
		//cnt = 0;
		psum_table.resize(transformed_.size());
		w.psum_table.resize(transformed_.size());
		csum_table.resize(transformed_.size());
		w.csum_table.resize(transformed_.size());

		for (int i = 0; i < psum_table.size(); ++i) {
			psum_table[i] = 1e255;
			w.psum_table[i] = 1e255;
			csum_table[i] = 1e255;
			w.csum_table[i] = 1e255;
		}

		std::vector<double> result(transformed_.size());
		// case 1
		result[0] += transformed_[0] * w.transformed_[0];
		
		int size = 1;
		// case 2

		// 普通に
		for (int level = 0; ; ++level) {
			const int Cuvw = 1 << level;
			for (int i = 0; i < size; ++i) {
				for (int j = 0; j < size; ++j) {
					int x10, y10;
					int x01, y01;
					int x11, y11;

					raw_index(level, i, j, 0x10, &x10, &y10);
					raw_index(level, i, j, 0x01, &x01, &y01);
					raw_index(level, i, j, 0x11, &x11, &y11);


					// 10と01の組み合わせ
					result[y11 * size_ + x11] += Cuvw * transformed_[y10 * size_ + x10] * w.transformed_[y01 * size_ + x01];
					result[y11 * size_ + x11] += Cuvw * transformed_[y01 * size_ + x01] * w.transformed_[y10 * size_ + x10];
					
					// 10と11の組み合わせ
					result[y01 * size_ + x01] += Cuvw * transformed_[y10 * size_ + x10] * w.transformed_[y11 * size_ + x11];
					result[y01 * size_ + x01] += Cuvw * transformed_[y11 * size_ + x11] * w.transformed_[y10 * size_ + x10];
					
					// 01と11の組み合わせ
					result[y10 * size_ + x10] += Cuvw * transformed_[y01 * size_ + x01] * w.transformed_[y11 * size_ + x11];
					result[y10 * size_ + x10] += Cuvw * transformed_[y11 * size_ + x11] * w.transformed_[y01 * size_ + x01];
				}
			}
			size *= 2;
			if (size == size_)
				break;
		}

		// case 3
		// ふつうに
		size = 1;
		result[0] += transformed_[1] * w.transformed_[1];
		result[0] += transformed_[size_] * w.transformed_[size_];
		result[0] += transformed_[size_ + 1] * w.transformed_[size_ + 1];
		result[0] += csum(0x00, w, *this, 0, 0, 0);
		for (int level = 0; ; ++level) {
			for (int i = 0; i < size; ++i) {
				for (int j = 0; j < size; ++j) {
					int x10, y10;
					int x01, y01;
					int x11, y11;

					raw_index(level, i, j, 0x10, &x10, &y10);
					raw_index(level, i, j, 0x01, &x01, &y01);
					raw_index(level, i, j, 0x11, &x11, &y11);

					result[y10 * size_ + x10] += transformed_[y10 * size_ + x10] * psum(w, level, i, j);
					result[y01 * size_ + x01] += transformed_[y01 * size_ + x01] * psum(w, level, i, j);
					result[y11 * size_ + x11] += transformed_[y11 * size_ + x11] * psum(w, level, i, j);

					result[y10 * size_ + x10] += w.transformed_[y10 * size_ + x10] * psum(*this, level, i, j);
					result[y01 * size_ + x01] += w.transformed_[y01 * size_ + x01] * psum(*this, level, i, j);
					result[y11 * size_ + x11] += w.transformed_[y11 * size_ + x11] * psum(*this, level, i, j);
					
					result[y10 * size_ + x10] += csum(0x10, w, *this, level, i, j);
					result[y01 * size_ + x01] += csum(0x01, w, *this, level, i, j);
					result[y11 * size_ + x11] += csum(0x11, w, *this, level, i, j);
				}
			}
			size *= 2;
			if (size == size_)
				break;
		}

		//std::cout << cnt << std::endl;
		transformed_ = result;
	}

	void print_transformed() {
		for (int y = 0; y < size_; y ++) {
			for (int x = 0; x < size_; x ++) {
				printf("%7.4f ", transformed_[y*size_ + x] );
				// std::cout << transformed_[y * size_ + x] << " ";
			}
			std::cout << std::endl;
		}
	}
	void print_reversed() {
		for (int y = 0; y < size_; y ++) {
			for (int x = 0; x < size_; x ++) {
				printf("%7.4f ", reversed_[y*size_ + x] );
				// std::cout << reversed_[y * size_ + x] << " ";
			}
			std::cout << std::endl;
		}
	}


	void reverse() {
		std::vector<double> output = transformed_;
		std::vector<double> output_tmp = transformed_;

		int now_size = 2;
		for (int i = 0; ; ++i) {
			if (!wavelet2Dreverse(&output[0], &output_tmp[0], now_size))
				break;
			now_size *= 2;
			output = output_tmp;
		}

		reversed_ = output;
	}

	void print_image_product(const Wavelet &w) {
		std::vector<double> tmp;
		for (int y = 0; y < size_; y ++) {
			for (int x = 0; x < size_; x ++) {
				tmp.push_back(image_[y * size_ + x] * w.image_[y * size_ + x]);
//				std::cout << image_[y * size_ + x] * w.image_[y * size_ + x] << " ";
			}
		}

		save("product.raw", tmp);
	}

	void print() {
		for (int y = 0; y < size_; y ++) {
			for (int x = 0; x < size_; x ++) {
				printf("%7.4f ", image_[y*size_ + x] );
				//;std::cout << image_[y * size_ + x] << " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	double diff_transformed(const Wavelet &w) {
		double d = 0;
		for (int i = 0; i < transformed_.size(); ++i)
			d += fabs(transformed_[i] - w.transformed_[i]);
		return d;
	}


	void create_mipmap(const Wavelet &w) {
		std::vector<std::vector<double> >mipmap;

		std::vector<double> product(image_.size());
		
		int lg = 0;
		for (int i = 0; i < image_.size(); ++i) {
			lg ++;
			product[i] = image_[i] * w.image_[i];
		}

		mipmap.push_back(product);
		int cnt = 1;
		for (int size = size_ / 2; size >= 1; size /= 2) {
			mipmap.push_back(std::vector<double>(size * size));

			for (int y = 0; y < size; ++y) {
				for (int x = 0; x < size; ++x) {
					lg ++;
					mipmap[cnt][y * size + x] = (
						mipmap[cnt - 1][(2 * y) * size + (2 * x)] +
						mipmap[cnt - 1][(2 * y + 1) * size + (2 * x)] +
						mipmap[cnt - 1][(2 * y) * size + (2 * x + 1)] +
						mipmap[cnt - 1][(2 * y + 1) * size + (2 * x + 1)]) * 0.25;
				}
			}
			cnt ++;
		}
	}
};

int main() {
	timeBeginPeriod(1);
	Wavelet w1("ibl3.raw");
	Wavelet w2("brdf.raw");
	Wavelet ans("ibl3xbrdf.raw");
	std::cout << std::fixed;
	/*
	w1.print_image_product(w2);
	return 0;
	*/
	/*
	w1.print();
	w2.print();
	ans.print();
	return 0;
	*/

	std::cout << "transform" << std::endl;
	w1.transform();
	w1.normalize();
	w1.cut(0.95);
	w1.sparsing();
	std::cout << "transform" << std::endl;
	w2.transform();
	w2.normalize();
	w2.cut(0.95);
	w2.sparsing();
	std::cout << "transform" << std::endl;
	ans.transform();
	ans.normalize();
	
	std::cout << "product" << std::endl;
	/*
	w1.print_transformed();
	std::cout << std::endl;
	w2.print_transformed();
	std::cout << std::endl;
	ans.print_transformed();
	std::cout << std::endl;
	std::cout << std::endl;
	*/

	//DWORD dw1 = timeGetTime();
	DWORD dw1 = timeGetTime();
	w1.product_sublinear(w2);
//	w1.product_linear(w2);
	DWORD dw2 = timeGetTime();
	std::cout << (dw2 - dw1) << " ms" << std::endl;
	//w1.print_transformed();
	
	
	DWORD dw3 = timeGetTime();
	w1.create_mipmap(w2);
	DWORD dw4 = timeGetTime();
	std::cout << (dw4 - dw3) << " ms" << std::endl;
	
	std::cout << "diff: " << ans.diff_transformed(w1) << std::endl;
	
	w1.unnormalize();
	w1.reverse();
	w1.save_reversed("result.raw");
	

	timeEndPeriod(1);
}