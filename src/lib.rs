mod img_funs;
mod building_blocks;
mod snd_funs;
mod eeg_19c_funs;

pub use crate::img_funs::*;
pub use crate::building_blocks::*;
pub use crate::snd_funs::*;
pub use crate::eeg_19c_funs::*;


#[cfg(test)]
mod tests {
use super::*;

    #[test]
    fn wtf() {assert!(true);}

    #[test]
    fn img_sorta() {
        use noise::Fbm;
        use noise::utils::*;
        use noise::MultiFractal;
        use itertools::Itertools;
        
        let mut dur = std::time::Duration::from_millis(0);
        for x in 0..10 {
            let clamp = |x, lower, upper| if x < lower {lower} else if x > upper {upper} else {x};
            let size_x = 1024;
            let size_y = size_x;
            let fbm = Fbm::new().set_octaves(32).set_persistence(x as f64 / 10.);
            let pmb = PlaneMapBuilder::new(&fbm).set_size(size_x as usize, size_y as usize).build();
            let pix_span =     (0..size_x as usize).cartesian_product(0..size_y as usize);
            let mut result = image::GrayImage::new(size_x as u32, size_y as u32);
            pix_span.clone()
                .map(|(xx, yy)| pmb.get_value(xx, yy))
                .map(|val| (val+1.0)*128.0)
                .map(|val| (clamp(val, 0.0, 255.99)) as u8)
                .zip(pix_span)
                .for_each(|(val, (xx, yy))| result.put_pixel(xx as u32, yy as u32, image::Luma([val])));
            drop((pmb, fbm));
            
            let begin = std::time::Instant::now();
            let (hrm, _, _) = measure_lum_2d_u8_image(result);
            dur += begin.elapsed();
            dbg!(hrm);
        }
        println!("Time elapsed: {:?}", dur);
        
        use noise::Seedable;
        let mut dur = std::time::Duration::from_millis(0);
        for x in 0..10 {
            let clamp = |x, lower, upper| if x < lower {lower} else if x > upper {upper} else {x};
            let size_x = 1<<10;
            let size_y = size_x;
            let fbm_1 = Fbm::new().set_octaves(32).set_persistence(x as f64 / 10.).set_seed(0);
            let fbm_2 = Fbm::new().set_octaves(32).set_persistence(x as f64 / 10.).set_seed(1);
            let fbm_3 = Fbm::new().set_octaves(32).set_persistence(x as f64 / 10.).set_seed(2);
            let pmb_1 = PlaneMapBuilder::new(&fbm_1).set_size(size_x as usize, size_y as usize).build();
            let pmb_2 = PlaneMapBuilder::new(&fbm_2).set_size(size_x as usize, size_y as usize).build();
            let pmb_3 = PlaneMapBuilder::new(&fbm_3).set_size(size_x as usize, size_y as usize).build();
            let pix_span =     (0..size_x as usize).cartesian_product(0..size_y as usize);
            let mut result = image::RgbImage::new(size_x as u32, size_y as u32);
            let map_fn = |val: f64| (val+1.0)*128.0;
            let map_fn = |val: f64| clamp(map_fn(val), 0.0, 255.99) as u8;
            let map_fn = |x: [f64; 3]| [map_fn(x[0]), map_fn(x[1]), map_fn(x[2])];
            pix_span.clone()
                .map(|(xx, yy)| [pmb_1.get_value(xx, yy), pmb_2.get_value(xx, yy), pmb_3.get_value(xx, yy)])
                .map(map_fn)
                .zip(pix_span)
                .for_each(|(val, (xx, yy))| result.put_pixel(xx as u32, yy as u32, image::Rgb(val)));
            drop((pmb_1, fbm_1, pmb_2, fbm_2, pmb_3, fbm_3));
            
            let begin = std::time::Instant::now();
            let (hrm, _, _) = measure_rgb_2d_u8_image(result);
            dur += begin.elapsed();
            dbg!(hrm);
        }
        println!("Time elapsed: {:?}", dur);
    } 
}   
