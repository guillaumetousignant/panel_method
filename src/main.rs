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

    let mut cl_alpha: Vec<CLALPHA> = Vec::with_capacity(filenames.len());

    let psi_params = PLOTPARAMS {
        res_x: 120,
        res_y: 80,
        origin_x: -1.0,
        origin_y: -1.0,
        span_x: 3.0,
        span_y: 2.0,
    };

    let uxuy_params = PLOTPARAMS {
        res_x: 1,
        res_y: 80,
        origin_x: 0.25,
        origin_y: -1.0,
        span_x: 0.0,
        span_y: 2.0,
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
        let psi = calculate_psi(sin_alpha, cos_alpha, &psi_params, &bod, &cof);
        let uxuy = velocity_slice(sin_alpha, cos_alpha, &uxuy_params, &bod, &cof);

        println!("Time elapsed: {}ms.", now.elapsed().as_millis());
        print(&bod, &cpd, cl, cm);
        write_cp(&bod, &cpd);
        write_psi(&bod, &psi);
        write_uxuy(&bod, &uxuy);

        cl_alpha.push(CLALPHA {
            alpha: bod.alpha,
            cl,
            cm,
        });
    }    
    write_cl_alpha(&cl_alpha);
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

struct CLALPHA {
    alpha: f64,
    cl: f64,
    cm: f64,
}

struct UXUY {
    x: f64,
    y: Vec<f64>,
    u_x: Vec<f64>,
    u_y: Vec<f64>,
}

struct PSI {
    x: Vec<f64>,
    y: Vec<f64>,
    psi: Vec<f64>,
}

struct PLOTPARAMS {
    res_x: usize,
    res_y: usize,
    origin_x: f64,
    origin_y: f64,
    span_x: f64,
    span_y: f64,
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
                    (0.5 * ((dxjp*dxjp + dyjp*dyjp)/(dxj*dxj + dyj*dyj)).ln(), 
                        (dyjp*dxj - dxjp*dyj).atan2(dxjp*dxj + dyjp*dyj))
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
        cp,
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

fn calculate_psi(sin_alpha: f64, cos_alpha: f64, params: &PLOTPARAMS, bod: &BOD, cof: &COF) -> PSI {
    let x_vec: Vec<f64> = (0..params.res_x).map(|x| params.origin_x + params.span_x * (x as f64)/(params.res_x as f64)).collect();
    let y_vec: Vec<f64> = (0..params.res_y).map(|y| params.origin_y + params.span_y * (y as f64)/(params.res_y as f64)).collect();
    
    let qs = &cof.b[0..bod.ndtot];
    let gamma = cof.b[cof.n-1];
    let mut psi_mat = vec![0.0; x_vec.len()*y_vec.len()];

    for (i, x) in x_vec.iter().enumerate() {
        for (j, y) in y_vec.iter().enumerate() {
            psi_mat[i*y_vec.len() + j] = cos_alpha*y - sin_alpha*x;
            for (k, q) in qs.iter().enumerate() {
                psi_mat[i*y_vec.len() + j] -= q * y / ((x - bod.x_mid[k]).powi(2) + (y - bod.y_mid[k]).powi(2))
                                            + gamma * 0.5/std::f64::consts::PI * ((x - bod.x_mid[k]).powi(2) + (y - bod.y_mid[k]).powi(2)).sqrt().ln(); // No type deduction for pow??
            }
        }
    }

    PSI {
        x: x_vec,
        y: y_vec,
        psi: psi_mat,
    }
}

fn velocity_slice(sin_alpha: f64, cos_alpha: f64, params: &PLOTPARAMS, bod: &BOD, cof: &COF) -> UXUY {
    let y_vec: Vec<f64> = (0..params.res_y).map(|y| params.origin_y + params.span_y * (y as f64)/(params.res_y as f64)).collect();

    let q = &cof.b[0..bod.ndtot];
    let gamma = cof.b[cof.n-1];
    let mut u_x = vec![cos_alpha; y_vec.len()];
    let mut u_y = vec![sin_alpha; y_vec.len()];

    for (i, y) in y_vec.iter().enumerate() {
        for j in 0..bod.ndtot {
            let (flog, ftan) = match j == i {
                true => (0., std::f64::consts::PI),
                false => {
                    let dxj = params.origin_x - bod.x[j];
                    let dxjp = params.origin_x - bod.x[j+1];
                    let dyj = y - bod.y[j];
                    let dyjp = y - bod.y[j+1];
                    (0.5 * ((dxjp*dxjp + dyjp*dyjp)/(dxj*dxj + dyj*dyj)).ln(), 
                        (dyjp*dxj - dxjp*dyj).atan2(dxjp*dxj + dyjp*dyj))
                },
            };

            let aa = 0.5/std::f64::consts::PI * ftan;
            let b = 0.5/std::f64::consts::PI * flog;
            u_x[i] -= b*q[j] - gamma*aa;
            u_y[i] -= b*q[j] - gamma*aa; // Not sure
        }
    }

    UXUY {
        x: params.origin_x,
        y: y_vec,
        u_x,
        u_y,
    }
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

fn write_cp(bod: &BOD, cpd: &CPD) {
    let mut lines: Vec<String> = Vec::with_capacity(bod.ndtot+1);

    lines.push(format!("TITLE = \"CP at alpha= {}\"
VARIABLES = \"X\", \"Y\", \"Cp\", \"Ue\"
ZONE T= \"Zone     1\",  I= {},  J= 1,  DATAPACKING = POINT", bod.alpha, bod.ndtot));

    for (i, x_mid) in bod.x_mid.iter().enumerate() {
        lines.push(format!("{xmid:>12.9} {ymid:>12.9} {cp:>12.9} {ue:>12.9}",
                        xmid = x_mid,
                        ymid = bod.y_mid[i],
                        cp = cpd.cp[i],
                        ue = cpd.ue[i]));
    }

    let mut file = fs::File::create(format!("cp-{:0width$}.dat", bod.alpha, width = 3)).expect("Error, unable to create cp output file.");
    writeln!(file, "{}", lines.join("\n")).expect("Error, unable to write to cp output file.");
}

fn write_psi(bod: &BOD, psi: &PSI) {
    let mut lines: Vec<String> = Vec::with_capacity(psi.psi.len()+1);

    lines.push(format!("TITLE = \"Psi at alpha= {}\"
VARIABLES = \"X\", \"Y\", \"Psi\"
ZONE T= \"Zone     1\",  I= {},  J= {},  DATAPACKING = POINT", bod.alpha, psi.x.len(), psi.y.len()));

    for (i, x) in psi.x.iter().enumerate() {
        for (j, y) in psi.y.iter().enumerate() {
            lines.push(format!("{x:>12.9} {y:>12.9} {psi:>12.9}",
                        x = x,
                        y = y,
                        psi = psi.psi[i*psi.y.len() + j]));
        }
    }

    let mut file = fs::File::create(format!("psi-{:0width$}.dat", bod.alpha, width = 3)).expect("Error, unable to create psi output file.");
    writeln!(file, "{}", lines.join("\n")).expect("Error, unable to write to psi output file.");
}

fn write_uxuy(bod: &BOD, uxuy: &UXUY) {
    let mut lines: Vec<String> = Vec::with_capacity(uxuy.y.len()+1);

    lines.push(format!("TITLE = \"Velocity at c/4 at alpha= {}\"
VARIABLES = \"X\", \"Y\", \"U_x\", \"U_y\"
ZONE T= \"Zone     1\",  I= {},  J= 1,  DATAPACKING = POINT", bod.alpha, uxuy.y.len()));

    for (i, y) in uxuy.y.iter().enumerate() {
        lines.push(format!("{x:>12.9} {y:>12.9} {u_x:>12.9} {u_y:>12.9}",
                    x = uxuy.x,
                    y = y,
                    u_x = uxuy.u_x[i],
                    u_y = uxuy.u_y[i]));
    }

    let mut file = fs::File::create(format!("uxuy-{:0width$}.dat", bod.alpha, width = 3)).expect("Error, unable to create uxuy output file.");
    writeln!(file, "{}", lines.join("\n")).expect("Error, unable to write to uxuy output file.");
}

fn write_cl_alpha(cl_alpha: &Vec<CLALPHA>) {
    let mut lines: Vec<String> = Vec::with_capacity(cl_alpha.len()+1);

    lines.push(format!("TITLE = \"CL versus alpha\"
VARIABLES = \"alpha\", \"CL\", \"CM\"
ZONE T= \"Zone     1\",  I= {},  J= 1,  DATAPACKING = POINT", cl_alpha.len()));

    for run in cl_alpha {
        lines.push(format!("{alpha:>12.9} {cl:>12.9} {cm:>12.9}",
                        alpha = run.alpha,
                        cl = run.cl,
                        cm = run.cm
                        ));
    }

    let mut file = fs::File::create("clalpha.dat").expect("Error, unable to create cl vs alpha output file.");
    writeln!(file, "{}", lines.join("\n")).expect("Error, unable to write to cl vs alpha output file.");
}