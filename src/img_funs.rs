use crate::building_blocks::*;

#[cfg(feature = "parallel")]
use rayon::{iter::*, slice::*};

#[cfg(not(feature = "parallel"))]
pub fn measure_rgb_2d_u8_image (input: image::RgbImage) -> (f64, Vec<f64>, Vec<f64>) {
    let width = input.width() as usize;
    let height = input.height() as usize;

    let get_key_from_sample = |(x, y, rgb): (u32, u32, &image::Rgb<u8>)| -> u64 {
        let normalise_as_u8 = |q, max, min| ((q - min) * 256 / max) as u8;
        let norm_x = |x| normalise_as_u8(x, input.width(), 0);
        let norm_y = |y| normalise_as_u8(y, input.height(), 0);
        let arr = [norm_x(x), norm_y(y), rgb[0], rgb[1], rgb[2]];
        morton_encoding::morton_encode_u8_5d(arr)
    };
    
    let clzs = get_clzs(input.enumerate_pixels(), get_key_from_sample);
    let (mut tmp, mut lacun) = get_results_from_clzs(clzs, 64);
    tmp.drain(0..24);
    lacun.drain(0..24);
    finalise_results(tmp, lacun, width*height, 8, 40)
}

#[cfg(not(feature = "parallel"))]
pub fn measure_lum_2d_u8_image (input: image::GrayImage) -> (f64, Vec<f64>, Vec<f64>) {
    let width = input.width() as usize;
    let height = input.height() as usize;

    let get_key_from_sample = |(x, y, lum): (u32, u32, &image::Luma<u8>)| -> u32 {
        let normalise_as_u8 = |q, max, min| ((q - min) * 256 / max) as u8;
        let norm_x = |x| normalise_as_u8(x, input.width(), 0);
        let norm_y = |y| normalise_as_u8(y, input.height(), 0);
        let arr = [norm_x(x), norm_y(y), lum[0]];
        morton_encoding::morton_encode_u8_3d(arr)
    };
    
    let clzs = get_clzs(input.enumerate_pixels(), get_key_from_sample);
    let (mut tmp, mut lacun) = get_results_from_clzs(clzs, 32);
    tmp.drain(0..8);
    lacun.drain(0..8);
    finalise_results(tmp, lacun, width*height, 8, 24)
}

#[cfg(feature = "parallel")]
pub fn measure_lum_2d_u8_image (input: image::GrayImage) -> (f64, Vec<f64>, Vec<f64>) {
// y as usize*self.width as usize + x as usize
    let width = input.width() as usize;
    let height = input.height() as usize;
    let buf = input.into_raw().into_par_iter();
    
    /*let x_range = (0..width).cycle().take(width * height);
    let y_range = (0..height).flat_map(|x| (x..=x).cycle().take(width));
    let buf = buf.zip(x_range).zip(y_range);
    
    
    let pixel_range = (0..height)
        .into_par_iter()
        .flat_map(|y| (0..width).into_par_iter().map(move |x| (x, y)))
        .map(move |(x, y)| (x, y, buf[y*width + x]));;
        
        
    
    let buf = buf.par_chunks(width)
        .enumerate()
        .flat_map(|(y,y_slice)|{
            y_slice.into_par_iter()
             .enumerate()
             .map(move |(x, lum)| (x, y, lum))
        });
    */
    let buf = buf.enumerate();
    
    // let get_key_from_sample = |((x, y), lum): ((usize, usize), u8)| -> u32 {
    let get_key_from_sample = |(comb, lum): (usize, u8)| -> u32 {
        let normalise_as_u8 = |q, max, min| ((q - min) * 256 / max) as u8;
        let norm_x = |x| normalise_as_u8(x, width, 0);
        let norm_y = |y| normalise_as_u8(y, height, 0);
        let x = comb % width;
        let y = comb / width;
        let x = norm_x(x);
        let y = norm_y(y);
        let arr = [x, y, lum];
        morton_encoding::morton_encode_u8_3d(arr)
    };
    
    let clzs = get_clzs_par(buf, get_key_from_sample);
    let (mut tmp, mut lacun) = get_results_from_clzs(clzs, 32);
    tmp.drain(0..8);
    lacun.drain(0..8);
    finalise_results(tmp, lacun, width*height, 8, 24)
}

 
#[cfg(feature = "parallel")]
pub fn measure_rgb_2d_u8_image (input: image::RgbImage) -> (f64, Vec<f64>, Vec<f64>) {
    
    
// y as usize*self.width as usize + x as usize
    let width = input.width() as usize;
    let height = input.height() as usize;
    let buf = input.into_raw();
    let buf = buf.par_chunks(3).enumerate();
    
    let get_key_from_sample = |(comb, rgb): (usize, &[u8])| -> u64 {
        let normalise_as_u8 = |q, max, min| ((q - min) * 256 / max) as u8;
        let norm_x = |x| normalise_as_u8(x, width, 0);
        let norm_y = |y| normalise_as_u8(y, height, 0);
        let x = norm_x(comb % width);
        let y = norm_y(comb / width);
        let arr = [x, y, rgb[0], rgb[1], rgb[2]];
        morton_encoding::morton_encode_u8_5d(arr)
    };
    
    let clzs = get_clzs_par(buf, get_key_from_sample);
    let (mut tmp, mut lacun) = get_results_from_clzs(clzs, 64);
    tmp.drain(0..24);
    lacun.drain(0..24);
    
    finalise_results(tmp, lacun, width*height, 8, 40)
}
