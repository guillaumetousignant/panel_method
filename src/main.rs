// MCG 5136 Assignment 2
// Guillaume Tousignant, 0300151859
// February 3rd, 2020

use std::error::Error;
use std::fs;
use std::io::Write;
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let filenames = match args.len() > 1 {
        true => args[1..].to_vec(),
        false => vec![get_input("Enter input file name (include extension name)")],
    };

    for filename in &filenames {
        let now = Instant::now();    
        println!("Input file is  {}", filename);    

        let bod = read_file(filename);
        let cos_alpha  = bod.alpha.to_radians().cos();
        let sin_alpha  = bod.alpha.to_radians().sin();

        let mut cof = coef(sin_alpha, cos_alpha, &bod);
        gauss(1, &mut cof);
        let cpd = vpdis(sin_alpha, cos_alpha, &bod, &cof);
        let (cl, cm) = clcm(sin_alpha, cos_alpha, &bod, &cpd);

        println!("Time elapsed: {}ms.", now.elapsed().as_millis());
        print(&bod, &cpd, cl, cm);
        plot(&bod, &cpd);
    }    
}

struct BOD {
    ndtot: usize,
    x: Vec<f64>,
    y: Vec<f64>,
    x_mid: Vec<f64>,
    y_mid: Vec<f64>,
    alpha: f64,
    costhe: Vec<f64>,
    sinthe: Vec<f64>,
}

struct COF {
    a: Vec<f64>,
    b: Vec<f64>,
    n: usize,
}

struct CPD {
    ue: Vec<f64>,
    cp: Vec<f64>,
}

fn read_file(filename: &String) -> BOD {
    let data = match fs::read_to_string(&filename) {
        Err(why) => panic!("Couldn't open {}: {}", filename, why.description()),
        Ok(data) => data,
    };

    let mut words = data.split_whitespace();
    let ndtot: usize = match words.next() {
        Some(item) => item.parse().unwrap(),
        None => panic!("Error, {} empty.", filename),
    };

    let mut x = vec![0.; ndtot+1];
    let mut y = vec![0.; ndtot+1];
    let mut x_mid = vec![0.; ndtot];
    let mut y_mid = vec![0.; ndtot];
    let mut costhe = vec![0.; ndtot];
    let mut sinthe = vec![0.; ndtot];

    for i in 0..ndtot+1 {
        x[i] = match words.next() {
            Some(item) => item.parse().unwrap(),
            None => panic!("Error, end of file reached at {} before x got to {} in {}.", i, ndtot, filename),
        };
    }
    for i in 0..ndtot+1 {
        y[i] = match words.next() {
            Some(item) => item.parse().unwrap(),
            None => panic!("Error, end of file reached at {} before y got to {} in {}.", i, ndtot, filename),
        };
    }

    for i in 0..ndtot {
        x_mid[i] = 0.5 * (x[i] + x[i+1]);
        y_mid[i] = 0.5 * (y[i] + y[i+1]);
        let dx: f64 = x[i+1] - x[i];
        let dy: f64 = y[i+1] - y[i]; // sqrt not found if the type isn't specified on those...
        let dist = (dx*dx + dy*dy).sqrt();

        sinthe[i] = dy/dist;
        costhe[i] = dx/dist;
    }

    let alpha: f64 = match words.next() {
        Some(item) => item.parse().unwrap(),
        None => panic!("Error, end of file reached before alpha in {}.", filename),
    };

    BOD {
        ndtot,
        x,
        y,
        x_mid,
        y_mid,
        alpha,
        costhe,
        sinthe,
        }
}

fn coef(sin_alpha: f64, cos_alpha: f64, bod: &BOD) -> COF {
    let n = bod.ndtot + 1;
    let mut a = vec![0.; n*n];
    let mut bv = vec![0.; n];
    for j in 0..n {
        a[(n-1)*n + j] = 0.;
    }
    for i in 0..bod.ndtot {
        a[i*n + n-1] = 0.;
        for j in 0..bod.ndtot {
            let (flog, ftan) = match j == i {
                true => (0., std::f64::consts::PI),
                false => {
                    let dxj = bod.x_mid[i] - bod.x[j];
                    let dxjp = bod.x_mid[i] - bod.x[j+1];
                    let dyj = bod.y_mid[i] - bod.y[j];
                    let dyjp = bod.y_mid[i] - bod.y[j+1];
                    (0.5 * ((dxjp*dxjp + dyjp*dyjp)/(dxj*dxj + dyj*dyj)).ln(), (dyjp*dxj - dxjp*dyj).atan2(dxjp*dxj + dyjp*dyj))
                },
            };

            let ctimtj = bod.costhe[i]*bod.costhe[j] + bod.sinthe[i]*bod.sinthe[j];
            let stimtj  = bod.sinthe[i]*bod.costhe[j] - bod.costhe[i]*bod.sinthe[j];
            a[i*n + j] = 0.5/std::f64::consts::PI * (ftan*ctimtj + flog*stimtj);
            let b = 0.5/std::f64::consts::PI * (flog*ctimtj - ftan*stimtj);
            a[i*n + n-1] += b;

            if (i < 1) || (i >= bod.ndtot-1) {
                a[(n-1)*n + j] -= b;
                a[(n-1)*n + n-1] += a[i*n + j];
            }   
        }
        bv[i] = bod.sinthe[i]*cos_alpha - bod.costhe[i]*sin_alpha;
    }
    bv[n-1] = -(bod.costhe[0] + bod.costhe[bod.ndtot-1])*cos_alpha - (bod.sinthe[0] + bod.sinthe[bod.ndtot-1])*sin_alpha;
    
    COF {
        a, 
        b: bv, 
        n
    }
}   

// Assumes a is of size n*n and b is of size n
fn gauss(m: usize, cof: &mut COF) {
    let n = cof.n;

    for k in 0..n-1 {
        let kp = k + 1;
        for i in kp..n {
            let r = cof.a[i*n + k]/cof.a[k*n + k];
            for j in kp..n {
                cof.a[i*n + j] -= r*cof.a[k*n + j];
            }
            for j in 0..m {
                cof.b[i*m + j] -= r*cof.b[k*m + j];
            }
        }
    }


    for k in 0..m {
        cof.b[(n-1)*m + k] /= cof.a[(n-1)*n + n-1];
        for i in (0..n-1).rev() {
            let ip = i + 1;
            for j in ip..n {
                cof.b[i*m + k] -= cof.a[i*n + j] * cof.b[j*m + k];
            }
            cof.b[i*m + k] /= cof.a[i*n + i];
        }
    }
}

fn vpdis(sin_alpha: f64, cos_alpha: f64, bod: &BOD, cof: &COF) -> CPD {
    let q = &cof.b[0..bod.ndtot];
    let gamma = cof.b[cof.n-1];
    let mut cp = vec![0.; bod.ndtot];
    let mut ue = vec![0.; bod.ndtot];

    for i in 0..bod.ndtot {
        let mut v_tan = cos_alpha*bod.costhe[i] + sin_alpha*bod.sinthe[i];
        for j in 0..bod.ndtot {
            let (flog, ftan) = match j == i {
                true => (0., std::f64::consts::PI),
                false => {
                    let dxj = bod.x_mid[i] - bod.x[j];
                    let dxjp = bod.x_mid[i] - bod.x[j+1];
                    let dyj = bod.y_mid[i] - bod.y[j];
                    let dyjp = bod.y_mid[i] - bod.y[j+1];
                    (0.5 * ((dxjp*dxjp + dyjp*dyjp)/(dxj*dxj + dyj*dyj)).ln(), (dyjp*dxj - dxjp*dyj).atan2(dxjp*dxj + dyjp*dyj))
                },
            };

            let ctimtj = bod.costhe[i]*bod.costhe[j] + bod.sinthe[i]*bod.sinthe[j];
            let stimtj  = bod.sinthe[i]*bod.costhe[j] - bod.costhe[i]*bod.sinthe[j];
            let aa = 0.5/std::f64::consts::PI * (ftan * ctimtj + flog * stimtj);
            let b = 0.5/std::f64::consts::PI * (flog * ctimtj - ftan * stimtj);
            v_tan -= b*q[j] - gamma*aa;
        }
        cp[i]   = 1. - v_tan*v_tan;
        ue[i]   = v_tan;
    }

    CPD {
        ue, 
        cp
    }
}

fn clcm(sin_alpha: f64, cos_alpha: f64, bod: &BOD, cpd: &CPD) -> (f64, f64) {
    let mut cfx = 0.;
    let mut cfy = 0.;
    let mut cm = 0.;

    for i in 0..bod.ndtot {
        let dx = bod.x[i+1] - bod.x[i];
        let dy = bod.y[i+1] - bod.y[i];
        cfx += cpd.cp[i]*dy;
        cfy -= cpd.cp[i]*dx;
        cm += cpd.cp[i]*(dx*bod.x_mid[i] + dy*bod.y_mid[i]);
    }

    let cl = cfy*cos_alpha - cfx*sin_alpha;
    
    (cl, cm)
}

fn get_input(prompt: &str) -> String{
    println!("{}",prompt);
    let mut input = String::new();
    match std::io::stdin().read_line(&mut input) {
        Ok(_goes_into_input_above) => {},
        Err(_no_updates_is_fine) => {},
    }
    input.trim().to_string()
}

fn print(bod: &BOD, cpd: &CPD, cl: f64, cm: f64) {
    println!("\n SOLUTION AT ALPHA = {:>10.5}\n", bod.alpha);
    println!("    j    X(j)      Y(j)      CP(j)      UE(j)\n");

    for i in 0..bod.ndtot {
        println!("{index:>5} {xmid:>9.5} {ymid:>9.5} {cp:>9.5} {ue:>9.5}",
                index = i,
                xmid = bod.x_mid[i],
                ymid = bod.y_mid[i],
                cp = cpd.cp[i],
                ue = cpd.ue[i]);
    }

    println!("\n    CL = {CL:>9.5}    CM = {CM:>9.5}\n\n",
            CL = cl,
            CM = cm);
}

fn plot(bod: &BOD, cpd: &CPD) {
    let mut strings: Vec<String> = Vec::with_capacity(bod.ndtot+1);

    strings.push(format!("TITLE = \"CP along airfoil\"
VARIABLES = \"X\", \"Y\", \"Cp\", \"Ue\"
ZONE T= \"Zone     1\",  I= {},  J= 1,  DATAPACKING = POINT", bod.ndtot));

    for (i, x_mid) in bod.x_mid.iter().enumerate() {
        strings.push(format!("{xmid:>9.5} {ymid:>9.5} {cp:>9.5} {ue:>9.5}",
                        xmid = x_mid,
                        ymid = bod.y_mid[i],
                        cp = cpd.cp[i],
                        ue = cpd.ue[i]));
    }

    let mut file = fs::File::create(format!("cp-{}.dat", bod.alpha)).expect("Error, unable to create output file.");
    writeln!(file, "{}", strings.join("\n")).expect("Error, unable to write to output file.");
}