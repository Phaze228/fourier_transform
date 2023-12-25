const std = @import("std");
const math = std.math;
const Complex = math.complex.Complex;

const SAMPLE_COUNT: usize = math.pow(usize, 2, 10);
const SAMPLE_COUNT_FLOAT: f32 = math.pow(f32, 2.0, 10);

const I: Complex(f32) = math.complex.sqrt(Complex(f32).init(-1.0, 0));
const PI = std.math.pi;
const PI_CPLX: Complex(f32) = Complex(f32).init(PI, 0);

pub fn main() !void {
    // Start by generating our "Samples" by using a sin equation evauluated at fractions of the FFT_SIZE
    // Example being, sin(2*PI*t) where t is index/fft_size
    var the_samples = genSamples();
    var out_samps: [SAMPLE_COUNT]Complex(f32) = std.mem.zeroes([SAMPLE_COUNT]Complex(f32));
    var out_samps2: [SAMPLE_COUNT]Complex(f32) = std.mem.zeroes([SAMPLE_COUNT]Complex(f32));
    std.debug.print("--FFT ITERATIVE--\n", .{});
    timeThis(&the_samples, &out_samps, fft);
    std.debug.print("--DFT--\n", .{});
    timeThis(&the_samples, &out_samps2, dft);
    std.debug.print("EQUAL:: -> {any} \n", .{std.meta.eql(out_samps, out_samps2)});
    // std.debug.print("{s: ^5} : {s: ^18} : {s: ^19} \n", .{ "INDEX", "OUTPUT1", "OUTPUT2" });
    //
    for (out_samps, 0..) |samp, i| {
        const out_sin1 = samp.im;
        const out_sin2 = out_samps2[i].im;
        const out_cos1 = samp.re;
        const out_cos2 = out_samps2[i].re;
        std.debug.print("{d:^5} : {d:^8.5} - {d:^8.5} : {d:^8.5} - {d:^8.5}\n", .{ i, out_sin1, out_cos1, out_sin2, out_cos2 });
    }
    std.debug.print("\n\n\n", .{});
    const sample_amount: usize = 50;
    var input_wave: [sample_amount]f32 = std.mem.zeroes([sample_amount]f32);
    var dft_wave: [sample_amount]Complex(f32) = std.mem.zeroes([sample_amount]Complex(f32));
    print_wave(&input_wave, &dft_wave, 3.17, sample_amount);

    // const r: f32 = @floor(@log2(SAMPLE_COUNT_FLOAT));
    // const r_usize: usize = @intFromFloat(r);
    // const r_rev: usize = math.shl(usize, @bitReverse(r_usize), r_usize);
    // std.debug.print("Bits: {b:.20}\n Reversed: {b:.20}\n", .{ r_usize, (r_rev) });
    // ft.samples_cplx = out_samps;
    // ft.print_cplx_samples();
}

pub fn timeThis(in: []f32, out: []Complex(f32), ft: *const fn (in: []f32, out: []Complex(f32), fft_size: usize) void) void {
    var time_start = std.time.milliTimestamp();
    ft(in, out, SAMPLE_COUNT);
    var time_end = std.time.milliTimestamp();
    std.debug.print("TIME ELAPSED: {d} milliseconds \n", .{time_end - time_start});
}

pub fn hamming_window(samps: []f32, sample_count: f32) void {
    for (samps, 0..) |f, i| {
        const t: f32 = @as(f32, @floatFromInt(i)) / sample_count - 1;
        const hamming: f32 = 0.5 - (0.5 * math.cos(2 * math.pi * t));
        samps[i] = f * hamming;
    }
}

pub fn print_wave(inputs: []f32, fft_out: []Complex(f32), f: f32, sample_count: usize) void {
    var i: usize = 0;
    var j: f32 = 0;
    const sc_float: f32 = @floatFromInt(sample_count);
    while (i < sample_count) : (i += 1) {
        const k: f32 = @floatFromInt(i);
        const t: f32 = k / sc_float;
        inputs[i] = math.sin(2 * PI * f * t);
        // std.debug.print("{d:.5}\n", .{inputs[i]});
    }
    hamming_window(inputs, sc_float);
    dft(inputs, fft_out, sample_count);

    i = 0;
    var mA: f32 = 0;
    while (i < sample_count) : (i += 1) {
        const a = @fabs(fft_out[i].magnitude());
        mA = @max(mA, a);
    }
    i = 0;
    while (i < sample_count) : (i += 1) {
        const a = @fabs(fft_out[i].magnitude());
        const t = a / mA;
        j = 0;
        while (j < t * sc_float) : (j += 1) {
            std.debug.print("*", .{});
        }
        std.debug.print("\n", .{});
    }
    i = 0;
    while (i < sample_count) : (i += 1) {
        j = 0;
        const t: f32 = (inputs[i] + 1) / 2;
        while (j < t * sc_float) : (j += 1) {
            std.debug.print(" ", .{});
        }
        std.debug.print("*\n", .{});
    }
}

pub fn reverse_bits(v: usize, k: usize) usize {
    var c: usize = v;
    var l: usize = 0;
    var i: usize = 0;
    while (i < k) : (i += 1) {
        l = (l << 1) + (c & 1);
        c >>= 1;
    }
    return l;
}

// gen samples for whatever mixture of wave equations you need but t is always the interval of index/fft_size
pub fn genSamples() [SAMPLE_COUNT]f32 {
    var samples: [SAMPLE_COUNT]f32 = undefined;
    for (samples, 0..) |_, i| {
        const k: f32 = @floatFromInt(i);
        const t: f32 = k / SAMPLE_COUNT_FLOAT;
        const wave_equation: f32 = math.sin(2 * PI * t);
        samples[i] = wave_equation;
    }
    return samples;
}

pub fn printSamples(samples: []f32) void {
    std.debug.print("\nSAMPLES\n", .{});
    for (samples, 0..) |s, i| {
        std.debug.print("{d:.3} : {d:.5} \n", .{ i, s });
    }
}

pub fn printSamplesComplex(samples: []Complex(f32)) void {
    std.debug.print("\n COMPLEX SAMPLES \n", .{});
    std.debug.print("{s: ^6}| {s: ^9}| {s: ^9}\n", .{ "LEVEL", "COS", "SIN" });
    std.debug.print("________________________\n", .{});
    for (samples, 0..) |sample, i| {
        std.debug.print("{d: ^6}| {d: ^9.4}| {d: ^9.4} \n", .{ i, sample.re, sample.im });
    }
}

pub fn dft(input_samples: []f32, output_frequencies: []Complex(f32), sample_count: usize) void {
    const smpl_cnt_float: f32 = @floatFromInt(sample_count);
    for (input_samples, 0..) |_, frequency| {
        const f: f32 = @floatFromInt(frequency);
        const freq: Complex(f32) = Complex(f32).init(f, 0);
        var i: usize = 0;
        while (i < sample_count) : (i += 1) {
            const k: f32 = @floatFromInt(i);
            const t: Complex(f32) = Complex(f32).init(k / smpl_cnt_float, 0);
            // EULERS e^ix = cosx +isinx
            const cmplx_samples: Complex(f32) = Complex(f32).init(input_samples[i], 0);
            const x = t.mul(freq).mul(I).mul(PI_CPLX).mul(Complex(f32).init(2, 0));
            output_frequencies[frequency] = output_frequencies[frequency].add(cmplx_samples.mul(math.complex.exp(x)));
        }
    }
}

pub fn fft_recurse(input_samples: []f32, step: usize, output_frequencies: []Complex(f32), sample_count: usize) void {
    if (sample_count == 1) {
        output_frequencies[0] = Complex(f32).init(input_samples[0], 0);
        return;
    }
    var is = input_samples;
    var of = output_frequencies;
    fft_recurse(is, step * 2, output_frequencies, sample_count / 2);
    is.ptr += step;
    of.ptr += sample_count / 2;
    fft_recurse(is, step * 2, of, sample_count / 2);
    var k: usize = 0;
    while (k < sample_count / 2) : (k += 1) {
        const len = sample_count;
        const t: Complex(f32) = Complex(f32).init(@as(f32, @floatFromInt(k)) / @as(f32, @floatFromInt(len)), 0);
        const x = t.mul(I).mul(PI_CPLX).mul(Complex(f32).init(-2, 0));
        const cmplx_samples: Complex(f32) = output_frequencies[k + sample_count / 2];
        const val = cmplx_samples.mul(math.complex.exp(x));
        const e = output_frequencies[k];
        output_frequencies[k] = e.add(val);
        output_frequencies[k + len / 2] = e.sub(val);
    }
}

pub fn fft(input_samples: []f32, output_frequencies: []Complex(f32), sample_count: usize) void {
    var step: usize = 1;
    var r: usize = @intFromFloat(@floor(@log2(@as(f32, @floatFromInt(sample_count)))));
    var l: usize = 0;
    var k: usize = 0;
    var p: usize = 0;
    var q: usize = 0;
    while (k < sample_count) : (k += 1) {
        l = reverse_bits(k, r);
        output_frequencies[l] = Complex(f32).init(input_samples[k], 0);
    }

    var twiddle: Complex(f32) = undefined;
    var twid: Complex(f32) = undefined;
    k = 0;
    while (k < r) : (k += 1) {
        l = 0;
        const step_cplx = Complex(f32).init(@floatFromInt(step), 0);
        const euler: Complex(f32) = math.complex.exp(PI_CPLX.div(step_cplx).mul(I));
        twiddle = euler;
        while (l < sample_count) : (l += 2 * step) {
            // twiddle.im = @floor(twiddle.im);
            // std.debug.print("...............................\n", .{});
            // std.debug.print("Index: {d} | Twiddle Factor: {d} + {d}\n", .{ l, twiddle.re, twiddle.im });
            // std.debug.print("_______________________________________\n", .{});
            twid = Complex(f32).init(1, 0);
            var n: usize = 0;
            while (n < step) : (n += 1) {
                p = l + n;
                q = p + step;
                output_frequencies[q] = output_frequencies[p].sub(twid.mul(output_frequencies[q]));
                output_frequencies[p] = (output_frequencies[p].mul(Complex(f32).init(2, 0))).sub(output_frequencies[q]);
                twid = twid.mul(twiddle);
            }
        }
        step <<= 1;
    }
}
