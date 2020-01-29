fn main() {
    let filename = get_input("Enter input file name (include extension name)");

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
    let mut a = vec![0.; n];
    let mut bv = vec![0.; n];
    for j in 0..n {
        a[n*n + j] = 0.;
    }
    for i in 0..n {
        a[i*n + n] = 0.;
        for j in 0..n {
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
            a[i*n + n] += b;

            if (i < 1) || (i >= n-1) {
                a[i*n + j] -= b;
                a[n*n + n] += a[i*n + j];
            }   
        }
        bv[i] = bod.sinthe[i]*cos_alpha - bod.costhe[i]*sin_alpha;
    }
    bv[n] = -(bod.costhe[0] + bod.costhe[bod.ndtot])*cos_alpha - (bod.sinthe[0] + bod.sinthe[bod.ndtot])*sin_alpha;
    
    COF {
        a, 
        b: bv, 
        n
    }
}   

// Assumes a is of size n*n and b is of size n
fn gauss(m: usize, cof: &mut COF) {
    let n = cof.b.len();

    for k in 0..n-1 {
        let kp = k + 1;
        for i in kp..n {
            let r = cof.a[i*n + k]/cof.a[k*n + k];
            for j in kp..n {
                cof.a[i*n + j] -= r*cof.a[k*n + j];
            }
            for j in 0..m {
                cof.b[i*n + j] -= r*cof.b[k*n + j];
            }
        }
    }

    for k in 0..m {
        cof.b[n*n + k] /= cof.a[n*n + n];
        for i in (0..n-1).rev() {
            let ip = i + 1;
            for j in ip..n {
                cof.b[i*n + k] -= cof.a[i*n + j] * cof.b[j*n + k];
            }
            cof.b[i*n + k] /= cof.a[i*n + i];
        }
    }

}

fn vpdis(sin_alpha: f64, cos_alpha: f64, bod: &BOD, cof: &COF) -> CPD {
    let q = &cof.b[0..bod.ndtot];
    let gamma = cof.b[cof.n-1];
    let mut cp = vec![0.; bod.ndtot];
    let mut ue = vec![0.; bod.ndtot];

    println!("    j    X(j)      Y(j)      CP(j)      UE(j)\n");
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
            v_tan -= b*q[j] + gamma*aa;
        }
        cp[i]   = 1. - v_tan*v_tan;
        ue[i]   = v_tan;
        println!("{index} {xmid} {ymid} {cp} {ue}",
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
    println!("\n\n    CL = {CL}    CM = {CM}",
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