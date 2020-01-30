use std::error::Error;
use std::fs;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let filename = match args.len() > 1 {
        true => args[1].clone(),
        false => get_input("Enter input file name (include extension name)"),
    };
    
    println!("Input file is  {}", filename);

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
    println!("\n\n SOLUTION AT ALPHA = {:>10.5}\n", alpha);

    let cos_alpha  = alpha.to_radians().cos();
    let sin_alpha  = alpha.to_radians().sin();

    let bod = BOD {
                ndtot,
                x,
                y,
                x_mid,
                y_mid,
                costhe,
                sinthe,
                };

    let mut cof = coef(sin_alpha, cos_alpha, &bod);
    gauss(1, &mut cof);
    let cpd = vpdis(sin_alpha, cos_alpha, &bod, &cof);
    clcm(sin_alpha, cos_alpha, &bod, &cpd);
}

struct BOD {
    ndtot: usize,
    x: Vec<f64>,
    y: Vec<f64>,
    x_mid: Vec<f64>,
    y_mid: Vec<f64>,
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

fn coef(sin_alpha: f64, cos_alpha: f64, bod: &BOD) -> COF{
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
    //println!("gamma: {}", gamma); GOOD

    println!("    j    X(j)      Y(j)      CP(j)      UE(j)\n");
    for i in 0..bod.ndtot {
        let mut v_tan = cos_alpha*bod.costhe[i] + sin_alpha*bod.sinthe[i];
        //println!("{}  {}", i, v_tan); // GOOD
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
            //println!("{},{}  {}", i, j, aa); // GOOD
            //println!("{},{}  {}", i, j, b); // GOOD
            v_tan -= b*q[j] - gamma*aa;
        }
        cp[i]   = 1. - v_tan*v_tan;
        ue[i]   = v_tan;
        println!("{index:>5} {xmid:>9.5} {ymid:>9.5} {cp:>9.5} {ue:>9.5}",
                index = i,
                xmid = bod.x_mid[i],
                ymid = bod.y_mid[i],
                cp = cp[i],
                ue = ue[i]);
    }

    CPD {
        ue, 
        cp
    }
}

fn clcm(sin_alpha: f64, cos_alpha: f64, bod: &BOD, cpd: &CPD) {
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
    println!("\n\n    CL = {CL:>9.5}    CM = {CM:>9.5}",
            CL = cl,
            CM = cm);
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