#include "led-matrix.addon.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <math.h>

#define BILLION  1000000000L;
#define MILLION  1000000000L;

inline double get_now_ms() {
	struct timespec t;
	if (clock_gettime(CLOCK_MONOTONIC_RAW, &t) < 0) {
		throw "Failed to get the current time.";
	}

	return (t.tv_sec * 1000) + (t.tv_nsec / 1000000);
}

using namespace rgb_matrix;

Napi::FunctionReference LedMatrixAddon::constructor;

Napi::Object LedMatrixAddon::Init(Napi::Env env, Napi::Object exports) {
	Napi::Function func = DefineClass(
	  env,
	  "LedMatrix",
	  { StaticMethod("defaultMatrixOptions", &LedMatrixAddon::default_matrix_options),
		StaticMethod("defaultRuntimeOptions", &LedMatrixAddon::default_runtime_options),
		InstanceMethod("afterSync", &LedMatrixAddon::after_sync),
		InstanceMethod("bgColor", &LedMatrixAddon::bg_color),
		InstanceMethod("brightness", &LedMatrixAddon::brightness),
		InstanceMethod("clear", &LedMatrixAddon::clear),
		InstanceMethod("drawBuffer", &LedMatrixAddon::draw_buffer),
		InstanceMethod("drawCircle", &LedMatrixAddon::draw_circle),
		InstanceMethod("drawFilledCircle", &LedMatrixAddon::draw_filled_circle),
		InstanceMethod("drawLine", &LedMatrixAddon::draw_line),
		InstanceMethod("drawRect", &LedMatrixAddon::draw_rect),
		InstanceMethod("drawFilledRect", &LedMatrixAddon::draw_filled_rect),
		InstanceMethod("drawPolygon", &LedMatrixAddon::draw_polygon),
		InstanceMethod("drawFilledPolygon", &LedMatrixAddon::draw_filled_polygon),
		InstanceMethod("drawText", &LedMatrixAddon::draw_text),
		InstanceMethod("fgColor", &LedMatrixAddon::fg_color),
		InstanceMethod("fill", &LedMatrixAddon::fill),
		InstanceMethod("font", &LedMatrixAddon::font),
		InstanceMethod("getAvailablePixelMappers", &LedMatrixAddon::get_available_pixel_mappers),
		InstanceMethod("height", &LedMatrixAddon::height),
		InstanceMethod("luminanceCorrect", &LedMatrixAddon::luminance_correct),
		InstanceMethod("map", &LedMatrixAddon::map),
		InstanceMethod("pwmBits", &LedMatrixAddon::pwm_bits),
		InstanceMethod("setPixel", &LedMatrixAddon::set_pixel),
		InstanceMethod("sync", &LedMatrixAddon::sync),
		InstanceMethod("width", &LedMatrixAddon::width) });

	constructor = Napi::Persistent(func);
	constructor.SuppressDestruct();
	exports.Set("LedMatrix", func);

	return exports;
}

/**
 * Process matrix & runtime options and initialize the internal RGBMatrix.
 */
LedMatrixAddon::LedMatrixAddon(const Napi::CallbackInfo& info)
  : Napi::ObjectWrap<LedMatrixAddon>(info)
  , after_sync_cb_(Napi::FunctionReference())
  , fg_color_(Color(0, 0, 0))
  , bg_color_(Color(0, 0, 0))
  , font_(nullptr)
  , font_name_("")
  , t_start_(get_now_ms())
  , t_sync_ms_(0)
  , t_dsync_ms_(0) {
	auto env = info.Env();

	if (!info[0].IsObject()) {
		throw Napi::Error::New(env, "Constructor expects its first parameter to be an object of matrix options!");
	}
	if (!info[1].IsObject()) {
		throw Napi::Error::New(env, "Constructor expects its first parameter to be an object of runtime options!");
	}
	auto matrixOpts  = create_matrix_options(env, info[0].As<Napi::Object>());
	auto runtimeOpts = create_runtime_options(env, info[1].As<Napi::Object>());

	this->matrix_ = CreateMatrixFromOptions(matrixOpts, runtimeOpts);
	this->canvas_ = this->matrix_->CreateFrameCanvas();

	if (this->matrix_ == NULL) { throw Napi::Error::New(env, "Failed to create matrix."); }
}

LedMatrixAddon::~LedMatrixAddon(void) {
	std::cerr << "Destroying matrix" << std::endl;
	delete matrix_;
}

Napi::Value LedMatrixAddon::sync(const Napi::CallbackInfo& info) {
	const char* data;
	size_t len;

	canvas_->Serialize(&data, &len);
	canvas_ = matrix_->SwapOnVSync(canvas_);
	if (!canvas_->Deserialize(data, len)) {
		throw Napi::Error::New(info.Env(), "Failed to sync canvas buffer with matrix.");
	}

	auto env = info.Env();

	if (!after_sync_cb_.IsEmpty()) {
        auto now = get_now_ms();
        auto now_ms = now - t_start_;
        t_dsync_ms_ = now_ms - t_sync_ms_;
        t_sync_ms_ = now_ms;

		auto resync = after_sync_cb_.Call(info.This(), {
			info.This(),
			Napi::Number::New(env, t_dsync_ms_),
			Napi::Number::New(env, t_sync_ms_)
		});

        if (resync.ToBoolean() == true) {
            sync(info);
        }
	}

	return Napi::Number::New(env, 0);
}

Napi::Value LedMatrixAddon::after_sync(const Napi::CallbackInfo& info) {
	auto cb = info[0].As<Napi::Function>();

	assert(cb.IsFunction());

	after_sync_cb_ = Napi::Persistent(cb);
	after_sync_cb_.SuppressDestruct();

	return info.This();
}

Napi::Value LedMatrixAddon::map(const Napi::CallbackInfo& info) {
	auto cb = info[0].As<Napi::Function>();

	assert(cb.IsFunction());

	auto env = info.Env();
	auto now = get_now_ms();
	auto now_ms = Napi::Number::New(env, now - t_start_);

    Napi::Array coord_array = Napi::Array::New(env, 3);
    uint32_t zero = 0; // The compiler can't match the overloaded signature if given 0 explicitly
    uint32_t one = 1;
    uint32_t two = 2;

    auto i = 0;

    for (int x = 0; x < this->matrix_->width(); x++) {
        coord_array.Set(zero, x);

        for (int y = 0; y < this->matrix_->height(); y++) {
            coord_array.Set(one, y);
            coord_array.Set(two, i++);

            auto color = cb.Call(info.This(), {
                coord_array,
                now_ms
            });

            assert(color.IsNumber());

            const auto hex = color.As<Napi::Number>().Uint32Value();
            this->matrix_->SetPixel(x, y, 0xFF & (hex >> 16), 0xFF & (hex >> 8), 0xFF & hex);
        }
    }

	return info.This();
}

Napi::Value LedMatrixAddon::brightness(const Napi::CallbackInfo& info) {
	if (info.Length() > 0 && info[0].IsNumber()) {
		auto brightness = info[0].As<Napi::Number>().Uint32Value();
		this->matrix_->SetBrightness(brightness);
		return info.This();
	}
	else {
		return Napi::Number::New(info.Env(), this->matrix_->brightness());
	}
}

Napi::Value LedMatrixAddon::clear(const Napi::CallbackInfo& info) {
	if (info.Length() > 0) {
		const auto x0	= info[0].As<Napi::Number>().Uint32Value();
		const auto y0	= info[1].As<Napi::Number>().Uint32Value();
		const auto x1	= info[2].As<Napi::Number>().Uint32Value();
		const auto y1	= info[3].As<Napi::Number>().Uint32Value();
		const auto black = Color(0, 0, 0);
		for (auto y = y0; y <= y1; y++) { DrawLine(this->canvas_, x0, y, x1, y, black); }
	}
	else {
		this->canvas_->Clear();
	}
	return info.This();
}

Napi::Value LedMatrixAddon::draw_buffer(const Napi::CallbackInfo& info) {
	const auto buffer = info[0].As<Napi::Buffer<uint8_t> >();
	const auto w	  = info[1].IsNumber() ? info[1].As<Napi::Number>().Uint32Value() : this->matrix_->width();
	const auto h	  = info[2].IsNumber() ? info[2].As<Napi::Number>().Uint32Value() : this->matrix_->height();
	const auto data   = buffer.Data();
	const auto len	= buffer.Length();

	assert(len == w * h * 3);

	Image* img	= new Image();
	Pixel* pixels = (Pixel*) malloc(sizeof(Pixel) * w * h);
	for (unsigned int i = 0; i < w * h; i++) {
		auto j = i * 3;
		Pixel p;
		p.r(data[j]);
		p.g(data[j + 1]);
		p.b(data[j + 2]);
		pixels[i] = p;
	}

	img->setPixels(w, h, pixels);

	assert(img->isValid());

	for (unsigned int y = 0; y < h; y++) {
		if (y > h) break;
		for (unsigned int x = 0; x < w; x++) {
			if (x > w) break;
			auto pixel = img->getPixel(x, y);
			this->canvas_->SetPixel(x, y, pixel.r(), pixel.g(), pixel.b());
		}
	}

	delete img;

	return info.This();
}

Napi::Value LedMatrixAddon::draw_circle(const Napi::CallbackInfo& info) {
	const auto x = info[0].As<Napi::Number>().Uint32Value();
	const auto y = info[1].As<Napi::Number>().Uint32Value();
	const auto r = info[2].As<Napi::Number>().Uint32Value();
	DrawCircle(this->canvas_, x, y, r, fg_color_);

	return info.This();
}

Napi::Value LedMatrixAddon::draw_filled_circle(const Napi::CallbackInfo& info) {
	const auto x = info[0].As<Napi::Number>().Int32Value();
	const auto y = info[1].As<Napi::Number>().Int32Value();
	const auto r = info[2].As<Napi::Number>().Int32Value();
	DrawCircle(this->canvas_, x, y, r, fg_color_);

	for(int i = 0 - r; i <= r; i++){
		for(int j = 0 -  r; j <= r; j++){
			if(pow(i, 2) + pow(j, 2) < pow(r, 2)){
				this->canvas_->SetPixel(i+x, y-j, fg_color_.r, fg_color_.g, fg_color_.b);
			}
		}
	}

	return info.This();
}

Napi::Value LedMatrixAddon::draw_line(const Napi::CallbackInfo& info) {
	const auto x0 = info[0].As<Napi::Number>().Uint32Value();
	const auto y0 = info[1].As<Napi::Number>().Uint32Value();
	const auto x1 = info[2].As<Napi::Number>().Uint32Value();
	const auto y1 = info[3].As<Napi::Number>().Uint32Value();
	DrawLine(this->canvas_, x0, y0, x1, y1, fg_color_);

	return info.This();
}

Napi::Value LedMatrixAddon::draw_rect(const Napi::CallbackInfo& info) {
	const auto x0 = info[0].As<Napi::Number>().Uint32Value();
	const auto y0 = info[1].As<Napi::Number>().Uint32Value();
	const auto w  = info[2].As<Napi::Number>().Uint32Value();
	const auto h  = info[3].As<Napi::Number>().Uint32Value();

	DrawLine(this->canvas_, x0, y0, x0 + w, y0, fg_color_);
	DrawLine(this->canvas_, x0 + w, y0, x0 + w, y0 + h, fg_color_);
	DrawLine(this->canvas_, x0 + w, y0 + h, x0, y0 + h, fg_color_);
	DrawLine(this->canvas_, x0, y0 + h, x0, y0, fg_color_);

	return info.This();
}

Napi::Value LedMatrixAddon::draw_filled_rect(const Napi::CallbackInfo& info) {
	const auto x0 = info[0].As<Napi::Number>().Int32Value();
	const auto y0 = info[1].As<Napi::Number>().Int32Value();
	const auto w  = info[2].As<Napi::Number>().Int32Value();
	const auto h  = info[3].As<Napi::Number>().Int32Value();

	for(int i = y0 + 1; i < y0 + h; i++){
		DrawLine(this->canvas_, x0+1, i, x0 + w -1, i, fg_color_);

	}

	DrawLine(this->canvas_, x0, y0, x0 + w - 1, y0, fg_color_); // Top
	DrawLine(this->canvas_, x0 + w - 1, y0, x0 + w - 1, y0 + h - 1, fg_color_); // Right
	DrawLine(this->canvas_, x0 + w - 1, y0 + h - 1, x0, y0 + h - 1, fg_color_);
	DrawLine(this->canvas_, x0, y0 + h - 1, x0, y0, fg_color_);

	return info.This();
}

Napi::Value LedMatrixAddon::draw_polygon(const Napi::CallbackInfo& info) {

	const auto coordinates = info[0].As<Napi::Array>();
	const auto length = coordinates.Length();
	
	for(int p = 0; p < length; p+=2){ // Iterate through points
		DrawLine(this->canvas_, coordinates[(p + 0) % length].As<Napi::Number>().Int32Value(), coordinates[(p + 1) % length].As<Napi::Number>().Int32Value(), coordinates[(p + 2) % length].As<Napi::Number>().Int32Value(), coordinates[(p + 3) % length].As<Napi::Number>().Int32Value(), fg_color_);
	}

	return info.This();
}

bool compareDoublesEqual(double a, double b){
	return std::abs(a - b) < (double)0.000001;
}

Napi::Value LedMatrixAddon::draw_filled_polygon(const Napi::CallbackInfo& info) {

	const auto coordinates = info[0].As<Napi::Array>();
	const auto length = coordinates.Length();
	std::vector<std::vector<double>> lines; // x0, y0, x1, y1, m, b

	int min_x = coordinates[(uint32_t)0].As<Napi::Number>().DoubleValue();
	int min_y = coordinates[(uint32_t)1].As<Napi::Number>().DoubleValue();
	int max_x = coordinates[(uint32_t)0].As<Napi::Number>().DoubleValue();
	int max_y = coordinates[(uint32_t)1].As<Napi::Number>().DoubleValue();

	int width;
	int height;

	double x0;
	double y0;
	double x1;
	double y1;
	double m;
	double b;

	bool fill_flag = false; // Running flag for filling polygon point-by-point.

	for(int p = 0; p < length; p+=2){ // Iterate through points and draw the lines
		x0 = (double)round(coordinates[(p + 0) % length].As<Napi::Number>().DoubleValue()); // Modulus helps us go back to the first point at the end of the loop.
		y0 = (double)round(coordinates[(p + 1) % length].As<Napi::Number>().DoubleValue());
		x1 = (double)round(coordinates[(p + 2) % length].As<Napi::Number>().DoubleValue());
		y1 = (double)round(coordinates[(p + 3) % length].As<Napi::Number>().DoubleValue());

		m = (y1 - y0) / (x1 - x0);
		b = y0 - (m * x0);

		lines.push_back(std::vector<double> {x0, y0, x1, y1, m, b});

		// Generate rectangle around polygon.
		min_x = fmin(fmin(min_x, x0), x1);
		min_y = fmin(fmin(min_y, y0), y1);
		max_x = fmax(fmax(max_x, x0), x1);
		max_y = fmax(fmax(max_y, y0), y1);

		DrawLine(this->canvas_, x0, y0, x1, y1, fg_color_);
	}

	width = (int)(max_x - min_x) + 1;
	height = (int)(max_y - min_y) + 1;

	for(int p = 0; p < width * height; p++){
		// Ray casting coordinates
		double y = min_y + floor(p / width);
		double x = min_x + (p % width);

		if(x == min_x){
			fill_flag = false;
		}

		std::vector<int> line_indexes_at_p; // Vector of index of line within "lines" vector.
		int lines_at_p = 0; // Like The Rentals' song.  :)

		// Loop through all lines at this point, p.
		for(int line = 0; line < lines.size(); line++){
			if(
				(
					(y == (int)round((lines[line][4] * x) + lines[line][5]) && lines[line][4] == 0) 
					|| 
					(x == (int)round((y - lines[line][5]) / lines[line][4]))
					|| isinf(lines[line][4])
				)
				&& (
					((x >= lines[line][0] && x <= lines[line][2]) || (x >= lines[line][2] && x <= lines[line][0])) 
					&& ((y >= lines[line][1] && y <= lines[line][3]) || (y >= lines[line][3] && y <= lines[line][1]))
				)
			){
				line_indexes_at_p.push_back(line);
				lines_at_p++;
			}
		}
		// Fill logic start

		// At a vertex vs. a line.
		if(lines_at_p > 1){
			// Here's where it gets nuts.
			// How many lines does the imaginary ray (y += 0.001) cross and if more than one, are the points on each line on same side of the imaginary line.
			// Calculate the imaginary x on each line of imaginary y.
			double imaginary_y = y + 0.001;
			double imaginary_x0 = (imaginary_y - lines[line_indexes_at_p[0]][5]) / lines[line_indexes_at_p[0]][4]; // Division by zero is happening.TODO
			double imaginary_x1 = (imaginary_y - lines[line_indexes_at_p[1]][5]) / lines[line_indexes_at_p[1]][4];

			// We have our imaginary y and our imaginary xes.  Are these valid points on each line?
			bool imaginary_0_valid_point = (
				(compareDoublesEqual(imaginary_y, (lines[line_indexes_at_p[0]][4] * imaginary_x0) + lines[line_indexes_at_p[0]][5]))
				&& (
					(imaginary_x0 >= lines[line_indexes_at_p[0]][0] && imaginary_x0 <= lines[line_indexes_at_p[0]][2]) 
					|| (imaginary_x0 <= lines[line_indexes_at_p[0]][0] && imaginary_x0 >= lines[line_indexes_at_p[0]][2])
				)
			);
			bool imaginary_1_valid_point = (
				(compareDoublesEqual(imaginary_y, (lines[line_indexes_at_p[1]][4] * imaginary_x1) + lines[line_indexes_at_p[1]][5]))
				&& (
					(imaginary_x1 >= lines[line_indexes_at_p[1]][0] && imaginary_x1 <= lines[line_indexes_at_p[1]][2]) 
					|| (imaginary_x1 <= lines[line_indexes_at_p[1]][0] && imaginary_x1 >= lines[line_indexes_at_p[1]][2])
				)
			);

			// If an imaginary point exists on each of the lines, toggle on and then off. (2 in my old thinking)
			if(imaginary_0_valid_point & imaginary_1_valid_point){
				this->canvas_->SetPixel(x, y, fg_color_.r, fg_color_.g, fg_color_.b);
			}
			// If an imaginary point exists on only one line (because of a cusp, for example), we've only crossed one line, effectively treating this as a side instead of a vertex.  Toggle fill on or off. (1 in my old thinking)
			if(imaginary_0_valid_point ^ imaginary_1_valid_point){
				fill_flag = !fill_flag;
			}
			// If an imaginary point exists on no line, do not toggle.  (0 in my old thinking)
			// if(!imaginary_0_valid_point & !imaginary_1_valid_point){
			// }
		}

		// At a line vs. a vertex.
		if(lines_at_p == 1 && lines[line_indexes_at_p[0]][4] != 0){ // And slope != 0
			fill_flag = !fill_flag;
		}

		if(fill_flag){
			this->canvas_->SetPixel(x, y, fg_color_.r, fg_color_.g, fg_color_.b);
		}
		// Fill logic end
	}

	return info.This();
}

Napi::Value LedMatrixAddon::draw_text(const Napi::CallbackInfo& info) {
	if (!font_) { throw Napi::Error::New(info.Env(), "Cannot draw text because the font has not been set!"); }
	const auto text		= std::string(info[0].As<Napi::String>()).c_str();
	const auto x		= info[1].As<Napi::Number>().Int32Value();
	const auto y		= info[2].As<Napi::Number>().Int32Value();
	const auto k		= info[3].IsNumber() ? info[3].As<Napi::Number>().Int32Value() : 0;
	const auto bg_color = bg_color_.r == 0 && bg_color_.g == 0 && bg_color_.b == 0 ? nullptr : &bg_color_;
	DrawText(this->canvas_, *font_, x, y + font_->baseline(), fg_color_, bg_color, text, k);

	return info.This();
}

Napi::Value LedMatrixAddon::fill(const Napi::CallbackInfo& info) {
	if (info.Length() > 0) {
		const auto x0 = info[0].As<Napi::Number>().Uint32Value();
		const auto y0 = info[1].As<Napi::Number>().Uint32Value();
		const auto x1 = info[2].As<Napi::Number>().Uint32Value();
		const auto y1 = info[3].As<Napi::Number>().Uint32Value();
		for (auto y = y0; y <= y1; y++) { DrawLine(this->canvas_, x0, y, x1, y, fg_color_); }
	}
	else {
		this->canvas_->Fill(fg_color_.r, fg_color_.g, fg_color_.b);
	}
	return info.This();
}

Napi::Value LedMatrixAddon::height(const Napi::CallbackInfo& info) {
	return Napi::Number::New(info.Env(), this->matrix_->height());
}

Napi::Value LedMatrixAddon::width(const Napi::CallbackInfo& info) {
	return Napi::Number::New(info.Env(), this->matrix_->width());
}

Napi::Value LedMatrixAddon::luminance_correct(const Napi::CallbackInfo& info) {
	if (info.Length() > 0 && info[0].IsBoolean()) {
		auto correct = info[0].As<Napi::Boolean>().ToBoolean();
		this->matrix_->set_luminance_correct(correct);
		return info.This();
	}
	else {
		return Napi::Boolean::New(info.Env(), this->matrix_->luminance_correct());
	}
}

Napi::Value LedMatrixAddon::pwm_bits(const Napi::CallbackInfo& info) {
	if (info.Length() > 0 && info[0].IsNumber()) {
		auto bits = info[0].As<Napi::Number>().Uint32Value();
		this->matrix_->SetPWMBits(bits);
		return info.This();
	}
	else {
		return Napi::Number::New(info.Env(), this->matrix_->pwmbits());
	}
}

Napi::Value LedMatrixAddon::set_pixel(const Napi::CallbackInfo& info) {
	const auto x = info[0].As<Napi::Number>().Uint32Value();
	const auto y = info[1].As<Napi::Number>().Uint32Value();
	this->canvas_->SetPixel(x, y, fg_color_.r, fg_color_.g, fg_color_.b);

	return info.This();
}

Napi::Value LedMatrixAddon::fg_color(const Napi::CallbackInfo& info) {
	if (info.Length() > 0) {
		auto color = LedMatrixAddon::color_from_callback_info(info);
		fg_color_  = color;
		return info.This();
	}
	else {
		return LedMatrixAddon::obj_from_color(info.Env(), fg_color_);
	}
}

Napi::Value LedMatrixAddon::bg_color(const Napi::CallbackInfo& info) {
	if (info.Length() > 0) {
		auto color = LedMatrixAddon::color_from_callback_info(info);
		bg_color_  = color;
		return info.This();
	}
	else {
		return LedMatrixAddon::obj_from_color(info.Env(), bg_color_);
	}
}

Napi::Value LedMatrixAddon::font(const Napi::CallbackInfo& info) {
	if (info.Length() > 0) {
		auto font   = Napi::ObjectWrap<FontAddon>::Unwrap(info[0].As<Napi::Object>());
		this->font_ = &(font->font);
		font_name_  = font->name(info).ToString();
		return info.This();
	}
	else {
		return Napi::String::New(info.Env(), font_name_);
	}
}

Napi::Value LedMatrixAddon::get_available_pixel_mappers(const Napi::CallbackInfo& info) {
    auto env = info.Env();
    auto mappers = GetAvailablePixelMappers();
    Napi::Array mapper_name_array = Napi::Array::New(env, mappers.size());

    for (uint8_t i = 0; i < mappers.size(); i++) {
        mapper_name_array.Set(i, Napi::String::New(env, mappers.at(i)));
    }

    return mapper_name_array;
}


/**
 * Create an instance of Options from a JS object.
 */
RGBMatrix::Options LedMatrixAddon::create_matrix_options(const Napi::Env& env, const Napi::Object& obj) {
	RGBMatrix::Options options = RGBMatrix::Options();

	options.brightness				 = obj.Get("brightness").As<Napi::Number>();
	options.chain_length			 = obj.Get("chainLength").As<Napi::Number>();
	options.cols					 = obj.Get("cols").As<Napi::Number>();
	options.disable_hardware_pulsing = obj.Get("disableHardwarePulsing").As<Napi::Boolean>();
	auto hardware_mapping			 = std::string(obj.Get("hardwareMapping").As<Napi::String>());
	options.hardware_mapping		 = strcpy(new char[hardware_mapping.size()], hardware_mapping.c_str());
	options.inverse_colors			 = obj.Get("inverseColors").As<Napi::Boolean>();
	auto led_rgb_sequence			 = std::string(obj.Get("ledRgbSequence").As<Napi::String>());
	options.led_rgb_sequence		 = strcpy(new char[led_rgb_sequence.size()], led_rgb_sequence.c_str());
	auto pixel_mapper_config		 = std::string(obj.Get("pixelMapperConfig").As<Napi::String>());
	options.pixel_mapper_config		 = strcpy(new char[pixel_mapper_config.size()], pixel_mapper_config.c_str());
	options.multiplexing			 = obj.Get("multiplexing").As<Napi::Number>();
	options.parallel				 = obj.Get("parallel").As<Napi::Number>();
	options.pwm_bits				 = obj.Get("pwmBits").As<Napi::Number>();
	options.pwm_dither_bits			 = obj.Get("pwmDitherBits").As<Napi::Number>();
	options.pwm_lsb_nanoseconds		 = obj.Get("pwmLsbNanoseconds").As<Napi::Number>();
	options.row_address_type		 = obj.Get("rowAddressType").As<Napi::Number>();
	options.rows					 = obj.Get("rows").As<Napi::Number>();
	options.scan_mode				 = obj.Get("scanMode").As<Napi::Number>();
	options.show_refresh_rate		 = obj.Get("showRefreshRate").As<Napi::Boolean>();

	// Validate the options using native method
	std::string error;
	if (!options.Validate(&error)) throw Napi::Error::New(env, error);
	return options;
}

/**
 * Create an instance of RuntimeOptions from a JS object.
 */
RuntimeOptions LedMatrixAddon::create_runtime_options(const Napi::Env& env, const Napi::Object& obj) {
	RuntimeOptions options = RuntimeOptions();

	options.gpio_slowdown   = obj.Get("gpioSlowdown").As<Napi::Number>();
	options.daemon			= obj.Get("daemon").As<Napi::Number>();
	options.drop_privileges = obj.Get("dropPrivileges").As<Napi::Number>();
	options.do_gpio_init	= obj.Get("doGpioInit").As<Napi::Boolean>();

	return options;
}

/**
 * Create a JS object from an instance of RGBMatrix::Options.
 */
Napi::Object LedMatrixAddon::matrix_options_to_obj(const Napi::Env& env, const RGBMatrix::Options& options) {
	auto obj = Napi::Object::New(env);

	std::string hardware_mapping = options.hardware_mapping == NULL ? "" : std::string(options.hardware_mapping);

	std::string led_rgb_sequence = options.led_rgb_sequence == NULL ? "" : std::string(options.led_rgb_sequence);

	std::string pixel_mapper_config
	  = options.pixel_mapper_config == NULL ? "" : std::string(options.pixel_mapper_config);

	obj.Set("brightness", Napi::Number::New(env, options.brightness));
	obj.Set("chainLength", Napi::Number::New(env, options.chain_length));
	obj.Set("cols", Napi::Number::New(env, options.cols));
	obj.Set("disableHardwarePulsing", Napi::Boolean::New(env, options.disable_hardware_pulsing));
	obj.Set("hardwareMapping", Napi::String::New(env, hardware_mapping));
	obj.Set("inverseColors", Napi::Boolean::New(env, options.inverse_colors));
	obj.Set("ledRgbSequence", Napi::String::New(env, led_rgb_sequence));
	obj.Set("multiplexing", Napi::Number::New(env, options.multiplexing));
	obj.Set("parallel", Napi::Number::New(env, options.parallel));
	obj.Set("pixelMapperConfig", Napi::String::New(env, pixel_mapper_config));
	obj.Set("pwmBits", Napi::Number::New(env, options.pwm_bits));
	obj.Set("pwmDitherBits", Napi::Number::New(env, options.pwm_dither_bits));
	obj.Set("pwmLsbNanoseconds", Napi::Number::New(env, options.pwm_lsb_nanoseconds));
	obj.Set("rowAddressType", Napi::Number::New(env, options.row_address_type));
	obj.Set("rows", Napi::Number::New(env, options.rows));
	obj.Set("scanMode", Napi::Number::New(env, options.scan_mode));
	obj.Set("showRefreshRate", Napi::Boolean::New(env, options.show_refresh_rate));

	return obj;
}

/**
 * Create a JS object from an instance of RuntimeOptions.
 */
Napi::Object LedMatrixAddon::runtime_options_to_obj(const Napi::Env& env, const RuntimeOptions& options) {
	auto obj = Napi::Object::New(env);

	obj.Set("gpioSlowdown", Napi::Number::New(env, options.gpio_slowdown));
	obj.Set("daemon", Napi::Number::New(env, options.daemon));
	obj.Set("dropPrivileges", Napi::Number::New(env, options.drop_privileges));
	obj.Set("doGpioInit", Napi::Boolean::New(env, options.do_gpio_init));

	return obj;
}

/**
 * Create a JS object from the default matrix options.
 */
Napi::Value LedMatrixAddon::default_matrix_options(const Napi::CallbackInfo& info) {
	auto env		   = info.Env();
	const auto options = RGBMatrix::Options();
	return LedMatrixAddon::matrix_options_to_obj(env, options);
}

/**
 * Create a JS object from the default runtime options.
 */
Napi::Value LedMatrixAddon::default_runtime_options(const Napi::CallbackInfo& info) {
	auto env = info.Env();
	return LedMatrixAddon::runtime_options_to_obj(env, RuntimeOptions());
}

/**
 * Create a Color instance from CallbackInfo.
 */
Color LedMatrixAddon::color_from_callback_info(const Napi::CallbackInfo& info) {
	if (info.Length() == 3) {
		uint8_t r = info[0].As<Napi::Number>().Uint32Value();
		uint8_t g = info[1].As<Napi::Number>().Uint32Value();
		uint8_t b = info[2].As<Napi::Number>().Uint32Value();
		return Color(r, g, b);
	}
	else if (info[0].IsObject()) {
		const auto obj = info[0].As<Napi::Object>();
		uint8_t r	  = obj.Get("r").As<Napi::Number>().Uint32Value();
		uint8_t g	  = obj.Get("g").As<Napi::Number>().Uint32Value();
		uint8_t b	  = obj.Get("b").As<Napi::Number>().Uint32Value();
		return Color(r, g, b);
	}
	else if (info[0].IsNumber()) {
		const auto hex = info[0].As<Napi::Number>().Uint32Value();
		return Color(0xFF & (hex >> 16), 0xFF & (hex >> 8), 0xFF & hex);
	}
	else {
		throw Napi::Error::New(info.Env(), "Failed to create color from parameters.");
	}
}

/**
 * Create an Object from a Color.
 */
Napi::Object LedMatrixAddon::obj_from_color(const Napi::Env& env, const Color& color) {
	Napi::Object obj = Napi::Object::New(env);
	obj.Set("r", color.r);
	obj.Set("g", color.g);
	obj.Set("b", color.b);
	return obj;
}
