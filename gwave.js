const earth_radius = 6370;
const am_dist_ratio = 1.07;
const pi=3.1415927;
const j = [0,1];
const speed_of_light = 299792.5;
const refractivity = .000315;
const radius_factor = 1.333333;

const c_in_air = speed_of_light / (1 + refractivity);
const adjusted_earth_radius = earth_radius*radius_factor;

//airy constants
const a0 = [2.3381074, 4.0879494, 5.5205598, 6.7867081, 7.9441336, 9.0226508, 10.0401743, 11.0085243, 11.9360156, 12.8287767];
const a1 = [1.0187930, 3.2481976, 4.8200992, 6.1633074, 7.3721773, 8.4884867, 9.5354491, 10.5276604, 11.4750566, 12.3847884];

//distances used to generate AM groundwave chart
const am_distances = require('./amDistances.json');

//conductivities used to generate AM groundwave chart
const am_sigmas = [0.1,0.5,1,1.5,2,3,4,5,6,7,8,10,15,20,30,40,5000];

//This maps assignable AM frequency ranges to their representative value for the charts.
const am_frequencies = require('./amFrequencies.json');

function realerfc(x) {
	var a1,a2,a3,a4,a5,p,t;
	a1 = .2548296;
	a2 = -.2844967;
	a3 = 1.4214137;
	a4 = -1.4531520;
	a5 = 1.0614054;
	p = .3275911;
	t = 1. / (1. + p*x);
	return t*(a1 + t*(a2 + t*(a3 + t*(a4 + t*(a5)))))*Math.exp(-Math.pow(x,2));
}

function salzer_f(x, y, n){
	return 2.*x - 2.*x*Math.cosh(n*y)*Math.cos(2.*x*y) + n*Math.sinh(n*y)*Math.sin(2.*x*y);
}

function salzer_g(x, y, n){
	return 2.*x*Math.cosh(n*y)*Math.sin(2.*x*y) + n*Math.sinh(n*y)*Math.cos(2.*x*y);
}

function salzer(z) {
	// double x = real(z);
	// double y = imag(z);
	const x = z[0];
	const y = z[1];


	var SUM=[0,0];
	if (y != 0) {
		const n = Math.min(Math.floor(Math.abs(80./y)), 50);
		// console.log(`Salzer x: ${x} y: ${y} n: ${n}`);
		for (var i=1; i<=n; i++) {
			const mult = Math.exp(-.25 * Math.pow(i,2)) / (Math.pow(i,2) + 4.*Math.pow(x,2));
			SUM = complexAdd(SUM, [mult * salzer_f(x,y,i), mult * salzer_g(x,y,i)]);  
			// console.log(`Salzer Sum ${i}: ${SUM}`);
		}
	}
	var salzer = complexMultiply(SUM, [-2. * Math.exp(-Math.pow(x,2)) / pi,0]);
	// console.log("Pre-Salzer:",salzer);
	if (x != 0) {
		// salzer -= Math.exp(-Math.pow(x,2)) * (1.-Math.cos(x*y) + j*Math.sin(x*y))/(pi*x);
		const sinxy = Math.sin(x*y);
		const cosxy = Math.cos(x*y);
		const w = [ //complex number
			Math.exp(-x*x) * sinxy / (pi * x),
			0
		];
		salzer = complexSubtract(salzer, complexMultiply(w, [sinxy, cosxy]));
	}else{
		// salzer -= j*y / pi;
		salzer = complexSubtract(salzer, [0, y/pi]);
	}

	if (x >= 0) {
		// salzer += realerfc(x);
		// console.log("realerfc(x): ", realerfc(x), "Salzer:", salzer);
		salzer = complexAdd(salzer, [realerfc(x),0]);
	}else{
		// salzer += 2. - realerfc(-x);
		// console.log("realerfc(-x): ", realerfc(-x), "Salzer:", salzer);
		salzer = complexAdd(salzer, [2-realerfc(-x), 0]);
	}
	return salzer;
}

function w_function(z) {
	const c1 = .4613135;
	const c2 = .09999216;
	const c3 = .002883894;
	const d1 = .1901635;
	const d2 = 1.7844927;
	const d3 = 5.5253437;
	var output = [0,0]; // complex
	output = complexAdd(output, w_factor(c1, z, d1));
	output = complexAdd(output, w_factor(c2, z, d2));
	output = complexAdd(output, w_factor(c3, z, d3));
	output = complexMultiply([0,1], complexMultiply(z,output));
	return output;
}

function w_factor(c, z, d) {
	return complexDivide([c,0], complexSubtract(complexMultiply(z,z), [d,0]));
}

function sommerfield(rho) {
	const rhoroot = complexSqrt(rho);
	// console.log(rhoroot);
	var za = [0,0];

	if (rhoroot[0] > 3.9 || rhoroot[1] > 3.0) {
		// 3rd method
		// console.log("3rd Sommerfield");
		var w = w_function(rhoroot);
		w = complexMultiply(w, rhoroot);
		w = complexMultiply([0,Math.sqrt(pi)], w);
		za = complexAdd([1,0], w);
	} else if (complexAbs(rhoroot) > 1) {
		// console.log("2nd Sommerfield");
		// 2nd method
		var erfc = salzer(complexMultiply([0,-1], rhoroot));
		erfc = complexMultiply(complexExp(complexMultiply([-1,0],rho)), erfc);
		erfc = complexMultiply(erfc, rhoroot);
		erfc = complexMultiply(erfc, [0,Math.sqrt(pi)]);
		za = complexAdd([1,0], erfc);
	} else {
		// 1st method
		var term = [1,0];
		var SUM = [1,0];
		for (var i=1; i<=33; i++) {
			// term = term * (-2 * rho) / (2 * i - 1);
			var A = [-2 / (2*i - 1),0];
			A = complexMultiply(rho, A);
			term = complexMultiply(term, A);
			// console.log(A, rho, term, SUM);
			SUM = complexAdd(SUM, term);
			if (complexAbs(term) < (complexAbs(SUM)/100000)) break;
		}
		za = complexAdd(SUM, complexMultiply(complexMultiply([0,1], rhoroot), complexMultiply(complexExp(complexMultiply([-1,0],rho)), [Math.sqrt(pi),0])));
	}
	return za;
}

function surface(rho, delta) {
	const za = sommerfield(rho);
	const rhoroot = complexSqrt(rho);
	var zadj;
	if (complexAbs(rho) > 0.5) {
		// console.log("Zadj3", za);
		// input.debug.surfaceMethod = "Zadj3";
		zadj = zadj3(delta, rho, za);
	}else{
		// console.log("Series");
		// input.debug.surfaceMethod = "Series";
		const z1 = complexMultiply([0,1],rhoroot);
		var z3 = complexPow(complexMultiply(delta,z1),3);
		// input.debug.z1 = z1;
		// input.debug.z3 = z3;
		// z3 = complexMultiply(complexMultiply(z3, z3), z3);
		zadj = complexSubtract(za, complexMultiply(z3,
			complexSubtract(series1(z1), complexMultiply(z3,
				complexSubtract(series2(z1), complexMultiply(z3, series3(z1)))
		))));
	}
	return complexAbs(zadj);
}

function zadj1(delta, rho, za) {
	const a = complexMultiply(za, complexAdd(complex(1),complexMultiply(complex(2),rho)));
	const b = complexSqrt(complexMultiply(complex(pi), rho));
	var c = complexSubtract(a,complex(1));
	// console.log(a,b,c);
	c = complexSubtract(c, complexMultiply(j, b));
	c = complexMultiply(c, complex(1/2));
	c = complexMultiply(c, complexPow(delta,3));
	return complexAdd(za, c);
}

function zadj2(delta, rho, za) {
	const a = complexMultiply(za, complexSubtract(complexDivide(complexPow(rho,2), complex(2)), complex(1)));
	var b = complexSubtract(complex(1), rho);
	b = complexMultiply(b, complexSqrt(complexMultiply(complex(pi), rho)));
	b = complexMultiply(j, b);
	var c = complexMultiply(complex(5/6), complexPow(rho,2));
	c = complexAdd(c, complex(1));
	c = complexSubtract(c, complexMultiply(complex(2), rho));
	var d = complexAdd(a,b);
	d = complexAdd(d,c);
	d = complexMultiply(d, complexPow(delta,6));
	return complexAdd(zadj1(delta,rho,za), d);
}

function zadj3(delta, rho, za) {
	var a = complexSubtract(complex(35/8), complexDivide(complexPow(rho,2), complex(4)));
	a = complexAdd(a, complexDivide(complexPow(rho,3), complex(6)));
	a = complexMultiply(za,a);
	var b = complexSubtract(complex(35/8), complexMultiply(complex(35/8),rho));
	b = complexAdd(b, complexMultiply(complex(31/16), complexPow(rho,2)));
	b = complexSubtract(b, complexMultiply(complex(5/16),complexPow(rho,3)));
	b = complexMultiply(b, complexMultiply(j, complexSqrt(complexMultiply(complex(pi),rho))));
	var c = complexSubtract(complex(35/8), complexMultiply(complex(35/4), rho));
	c = complexAdd(c, complexMultiply(complex(67/12),complexPow(rho,2)));
	c = complexSubtract(c, complexMultiply(complex(5/3),complexPow(rho,3)));
	var d = complexSubtract(a,b);
	d = complexSubtract(d,c);
	d = complexMultiply(d, complexPow(delta,9));
	return complexAdd(zadj2(delta,rho,za), d);
}

function series1(z) {
	const sqrtpi = Math.sqrt(3.1415927);
	var oddterm = complexMultiply(complex(4/(3*sqrtpi)),z);
	var eventerm = complex(1);
	var SUM = complexAdd(eventerm, complexMultiply(complex(2),oddterm));
	var term = complex(0);
	var TEST = complex(0);
	const zsquared = complexMultiply(z,z);
	//for i in range(2,50):
	for (var i=2; i<=50; i++) {
		if (i%2 == 0){
			// eventerm = 2.*eventerm*zsquared/(i + 2.);
			eventerm = complexMultiply(complex(2/(i+2)),eventerm);
			eventerm = complexMultiply(zsquared,eventerm);
			term = eventerm;
		}else{
			// oddterm = 2.*oddterm*zsquared/(i + 2.);
			oddterm = complexMultiply(complex(2/(i+2)),oddterm);
			oddterm = complexMultiply(zsquared,oddterm);
			term = oddterm;
		}
		// SUM += term * (i+1.);
		// console.log(i, eventerm, oddterm, SUM);
		SUM = complexAdd(SUM, complexMultiply(term,complex(i+1)));
		if (complexAbs(complexSubtract(SUM,TEST)) < (complexAbs(TEST)/1000000)) {
			break;
		}
		TEST=SUM;
	}
	return complexMultiply(SUM,complex(sqrtpi/2));
}

function series2(z) {
	const sqrtpi = Math.sqrt(3.1415927);
	var oddterm = complexDivide(z,complex(6));
	var eventerm = complexDivide(complex(8), complex(15*sqrtpi));
	var SUM = complexAdd(complexMultiply(complex(7),eventerm), complexMultiply(complex(2*8),oddterm));
	var term = complex(0);
	var TEST = complex(0);
	const zsquared = complexMultiply(z,z);
	//for i in range(2,50):
	for (var i=2; i<=50; i++) {
		if (i%2 == 0){
			// eventerm = 2.*eventerm*zsquared/(i + 2.);
			eventerm = complexMultiply(complex(2/(i+5)),eventerm);
			eventerm = complexMultiply(zsquared,eventerm);
			term = eventerm;
		}else{
			// oddterm = 2.*oddterm*zsquared/(i + 2.);
			oddterm = complexMultiply(complex(2/(i+5)),oddterm);
			oddterm = complexMultiply(zsquared,oddterm);
			term = oddterm;
		}
		// SUM += term * (i+1.);
		// console.log(i, eventerm, oddterm, SUM);
		SUM = complexAdd(SUM, complexMultiply(term,complex((i+1)*(i+7))));
		if (complexAbs(complexSubtract(SUM,TEST)) < (complexAbs(TEST)/1000000)) {
			break;
		}
		TEST=SUM;
	}
	return complexMultiply(SUM,complex(sqrtpi/8));
}

function series3(z) {
	const sqrtpi = Math.sqrt(3.1415927);
	var oddterm = complexMultiply(complex(32/(945*sqrtpi)),z);
	var eventerm = complex(1/24);
	var SUM = complexAdd(complexMultiply(complex(126),eventerm), complexMultiply(complex(2*147),oddterm));
	var term = complex(0);
	var TEST = complex(0);
	const zsquared = complexMultiply(z,z);
	//for i in range(2,50):
	for (var i=2; i<=50; i++) {
		if (i%2 == 0){
			// eventerm = 2.*eventerm*zsquared/(i + 2.);
			eventerm = complexMultiply(complex(2/(i+8)),eventerm);
			eventerm = complexMultiply(zsquared,eventerm);
			term = eventerm;
		}else{
			// oddterm = 2.*oddterm*zsquared/(i + 2.);
			oddterm = complexMultiply(complex(2/(i+8)),oddterm);
			oddterm = complexMultiply(zsquared,oddterm);
			term = oddterm;
		}
		// SUM += term * (i+1.);
		// console.log(i, eventerm, oddterm, SUM);
		SUM = complexAdd(SUM, complexMultiply(term,complex((i+1)*(i*i + 20*i+126))));
		if (complexAbs(complexSubtract(SUM,TEST)) < (complexAbs(TEST)/1000000)) {
			break;
		}
		TEST=SUM;
	}
	return complexMultiply(SUM,complex(sqrtpi/48));
}

function gwave(input) {

	input.debug = {};

	if (input.frequency > 2) input.frequency /= 1000;
	if (input.frequency > 2 || input.frequency < 0.5) return -1; 

	const wavelength = c_in_air / (input.frequency * 1e6);
	const x = 2*Math.pow(c_in_air*1e-5,2) * input.sigma/input.frequency;
	const b2 = Math.fround(Math.atan2(input.epsilon,x));
	const b1 = Math.fround(Math.atan2(input.epsilon-1,x));
	const b = 2*b2 - b1; //numerical distance phase
	const p = Math.fround(pi * input.distance * Math.pow(Math.cos(b2),2) / (wavelength * x * Math.cos(b1))); //numerical distance magnitude
	const k = Math.fround(Math.pow(wavelength / (2 * pi * adjusted_earth_radius), 1/3) * Math.sqrt(x * Math.cos(b1))/Math.cos(b2)); //Norton's K
	const rho = complexPolar(p,b);
	const rhoroot = complexSqrt(rho);
	const delta = complexPolar(k, pi * 3/4 - b/2);
	const chi = (input.distance / adjusted_earth_radius) * Math.pow(2*pi*adjusted_earth_radius/wavelength, 1/3);
	const critical_distance = 80 * Math.pow(1.3333, .6666) / Math.pow(input.frequency, .3333);

	var attenuation;

	if (input.distance < critical_distance) {
		attenuation = surface(rho, delta);
		// console.log("Surface Attenuation", attenuation);
		input.method = "Surface";
	}else{
		// console.log("Residues",input.distance, critical_distance);
		attenuation = residues(chi, k, b/2, delta);
		input.method = "Residues";
	}

	// console.log(`x: ${x}, b2: ${b2}, b1: ${b1}, b: ${b}, p: ${p}, k: ${k}, rho: ${rho}, rhoroot: ${rhoroot}, delta: ${delta}, chi: ${chi}, critical_distance: ${critical_distance}, attenuation: ${attenuation}`);

	input.attenuation = attenuation;
	input.rho = rho;
	input.rhoroot = rhoroot;
	input.delta = delta;
	input.x = x;
	input.b2 = b2;
	input.b1 = b1;
	input.b = b;
	input.p = p;
	input.k = k;
	input.chi = chi;
	input.critical_distance = critical_distance;
	input.output = (Math.sqrt(input.erp) * attenuation * input.field / input.distance);
	input.single_precision = Math.fround(input.output);
	return input;
}

function tfn(airy) {
// 	Function producing TAU0 or TAU1 from the amplitude AIRY0 or AIRY1:
	return complexMultiply(complex(airy / Math.pow(2,1/3)), complexExp([0,pi/3]));
}


function taufn0(tau, delta) {
	const c3 = complexMultiply(complex(-2/3),tau);
	const c5 = complexMultiply(complex(-4/5), complexPow(tau,2));
	const c6 = complexMultiply(complex(14/9),tau);
	// c7 = -(5.+8.* pow(tau,3.))/7.
	const c7 = complexDivide( complexAdd(complex(5), complexMultiply(complex(8), complexPow(tau,3))),complex(-7));
	// c8 = 58./15.* pow(tau,2.)
	const c8 = complexMultiply(complex(58/15), complexPow(tau,2));
	// c9 = -tau*(2296./567. +16./9. * pow(tau,3.) )
	var c9 = complexMultiply(complex(16/9), complexPow(tau,3));
	c9 = complexAdd(complex(2296/567), c9);
	c9 = complexMultiply(tau,c9);
	c9 = complexMultiply(complex(-1), c9);
	// c10 = 47./35. + 4656./525.* pow(tau,3.)
	var c10 = complexMultiply(complex(4656/525), complexPow(tau,3));
	c10 = complexAdd(complex(47/35), c10);

	var output = (c10);
	output = complexAdd(c9, complexMultiply(delta,output));
	output = complexAdd(c8, complexMultiply(delta,output));
	output = complexAdd(c7, complexMultiply(delta,output));
	output = complexAdd(c6, complexMultiply(delta,output));
	output = complexAdd(c5, complexMultiply(delta,output));
	output = complexAdd(complex(1/2), complexMultiply(delta,output));
	output = complexAdd(c3, complexMultiply(delta,output));
	output = complexAdd(complex(0), complexMultiply(delta,output));
	output = complexAdd(complex(-1), complexMultiply(delta,output));
	output = complexAdd(tau, complexMultiply(delta,output));
	return output;
}

function taufn1(tau, q) {
	//print(t,q)
	// complex<double> d1 = -1./(2.*t);
	const d1 = complexDivide(complex(-1/2), tau);
	// complex<double> d2 = -1./(8.*pow(t,3));
	const d2 = complexDivide(complex(-1/8),complexPow(tau,3));
	// complex<double> d3 = -(1./pow(t,2))*(1./12. + 1./(16.*pow(t,3)));
	var d3 = complexDivide(complex(1/16), complexPow(tau,3));
	d3 = complexAdd(complex(1/12),d3);
	d3 = complexMultiply(complexDivide(complex(-1),complexPow(tau,2)),d3);
	// complex<double> d4 = -(1.0/pow(t,4))*( 7.0/96.0 + 5.0/(128.0*pow(t,3) ) );
	var d4 = complexDivide(complex(5/128),complexPow(tau,3));
	d4 = complexAdd(complex(7/96), d4);
	d4 = complexMultiply(complexDivide(complex(-1),complexPow(tau,4)),d4);
	// complex<double> d5 = -(1.0/pow(t,3))*(1.0/40.0 + (1.0/pow(t,3))*(21.0/320. + 7.0/(256.0*pow(t,3)) ));
	var d5 = complexDivide(complex(7/56),complexPow(tau,3));
	d5 = complexAdd(complex(21/320), d5);
	d5 = complexMultiply(complexDivide(complex(1),complexPow(tau,3)),d5);
	d5 = complexAdd(complex(1/40), d5);
	d5 = complexMultiply(complexDivide(complex(-1),complexPow(tau,3)),d5);
	// complex<double> d6 = -(1.0/pow(t,5)) * ( 29.0/ 720.0 + 1.0/(pow(t,3)) * ( 77.0/1280.0 + 21.0/(1024.0*pow(t,3)) ));
	var d6 = complexDivide(complex(21/1024),complexPow(tau,3));
	d6 = complexAdd(complex(77/1280), d6);
	d6 = complexMultiply(complexDivide(complex(1),complexPow(tau,3)),d6);
	d6 = complexAdd(complex(29/720), d6);
	d6 = complexMultiply(complexDivide(complex(-1),complexPow(tau,5)),d6);
	// complex<double> d7 = -(1.0/pow(t,4))*(1.0/112.0 + (1.0/pow(t,3))*(19.0/360.0 + (1.0/pow(t,3))*( 143.0/2560 + 33.0/( 2048.0 * pow(t,3) ) )));
	var d7 = complexDivide(complex(33/2048),complexPow(tau,3));
	d7 = complexAdd(complex(143/2560), d7);
	d7 = complexMultiply(complexDivide(complex(1),complexPow(tau,3)),d7);
	d7 = complexAdd(complex(19/360), d7);
	d7 = complexMultiply(complexDivide(complex(1),complexPow(tau,3)),d7);
	d7 = complexAdd(complex(1/112), d7);
	d7 = complexMultiply(complexDivide(complex(-1),complexPow(tau,4)),d7);
	// complex<double> d8 = -(1.0/pow(t,6))*( 97.0/4480 + (1.0/pow(t,3))*(163.0/2560 + (1.0/pow(t,3))*(429.0/8192 + 429.0/( 32768.0 * pow(t,3) ) )));
	var d8 = complexDivide(complex(429/32768),complexPow(tau,3));
	d8 = complexAdd(complex(49/8192), d8);
	d8 = complexMultiply(complexDivide(complex(1),complexPow(tau,3)),d8);
	d8 = complexAdd(complex(163/2560), d8);
	d8 = complexMultiply(complexDivide(complex(1),complexPow(tau,3)),d8);
	d8 = complexAdd(complex(97/4480), d8);
	d8 = complexMultiply(complexDivide(complex(-1),complexPow(tau,6)),d8);
	// complex<double> out = t + q*(d1 + q*(d2 + q*(d3 + q*(d4 + q*(d5 + q*(d6 + q*(d7 + q*(d8))))))));

	var out = d8;
	out = complexAdd(d7, complexMultiply(q,out));
	out = complexAdd(d6, complexMultiply(q,out));
	out = complexAdd(d5, complexMultiply(q,out));
	out = complexAdd(d4, complexMultiply(q,out));
	out = complexAdd(d3, complexMultiply(q,out));
	out = complexAdd(d2, complexMultiply(q,out));
	out = complexAdd(d1, complexMultiply(q,out));
	out = complexAdd(tau, complexMultiply(q,out));

	return out;
}



function airy0(s) {
	 if(s < 10) {
	 	return a0[s];
	 }else{
		const x_s = 3. * pi * (4*(s+1) - 1)/8;
		return (Math.pow((x_s),2 / 3.0) * (1 + 5 / 48 * Math.pow(1.0 / (x_s),2.0)));
	}
}

function airy1(s) {
	if(s < 10) {
		return a1[s];
	}else{
		const y_s = (3.0 * pi * (4 * (s+1) - 3) / 8.0);
		return (Math.pow((y_s),2 / 3) * (1 - 7 / 48 * Math.pow(1.0 / (y_s),2.0)));
	}
}

function deltau(T, D, DELDEL) {
	var output = complexMultiply(T,complexPow(D,2));
	output = complexMultiply(complex(2), output);
	output = complexSubtract(output, complex(1));
	output = complexDivide(DELDEL, output);
	return output;
}

const residue_options = {
	precision: 1e-5,
	fineness: 1e-3,
	maxterms: 50,
	minsteps: 5
}

function residues(chi, k, psi, delta) {
	var deldel, delnew; //complex
	var test = complex(0);
	var zs = complex(0);
	var numpoints = 0;
	var tau=[];
	var qsqr = complexPow(complexDivide(complex(1),delta), 2);
	const N = Math.max(5, Math.floor(complexAbs(delta) / residue_options.fineness));
	// console.log(`Residues N: ${N}`);

	for (var s=1; s <= residue_options.maxterms; s++) {
		const tau0 = tfn(airy0(s-1));
		const tau1 = tfn(airy1(s-1));
		// console.log(`\tTau0: ${tau0}\tTau1: ${tau1}`);

		if ( complexAbs(complexMultiply(tau0,complex(Math.pow(k,2)))) < 0.16) {
			// console.log("Small k");
			tau[s] = taufn0(tau0, delta);
		}else if ( complexAbs(complexMultiply(tau1,complex(Math.pow(k,2)))) > 1.44) {
			// console.log("Large K");
			tau[s] = taufn1(tau1, complexDivide(complex(1), delta));
		}else{
			// console.log("Integrating");
			var T = tau0;
			// deldel = complexDivide(delta, complex(N));
			var z1;
			if ( complexAbs(delta) >= residue_options.fineness * residue_options.minsteps) {
				z1 = complexMultiply(complex(residue_options.fineness / complexAbs(delta)), delta);
			}else{
				z1 = complexDivide(delta, complex(residue_options.minsteps));
			}
			var del1 = complex(0);

			while (complexAbs(del1) - complexAbs(delta) < 0) {
				deldel = complexMultiply(z1, complex(Math.min(1,complexAbs(complexSubtract( complexMultiply(complex(2), complexMultiply(T, complexPow(del1,2))) , complex(1))))));
				T = taustep(T,del1,deldel);
				// console.log(T, del1, deldel);
				del1 = complexAdd(del1,deldel);
			}

			tau[s] = taustep(T, del1, complexSubtract(delta, del1));

			numpoints = s;
			// for (var i=1; i <= N+1; i++) {
			// 	const TK1 = deltau(T, del1, deldel);
			// 	del1 = complexAdd(del1, complexDivide(deldel,complex(2)));
			// 	const TK2 = deltau(complexAdd(T,complexDivide(TK1,complex(2))), del1, deldel);
			// 	const TK3 = deltau(complexAdd(T,complexDivide(TK2,complex(2))), del1, deldel);
			// 	del1 = complexAdd(del1, complexDivide(deldel,complex(2)));
			// 	const TK4 = deltau(complexAdd(T,TK3), del1, deldel);
			// 	T = complexAdd(T, complexDivide(complexAdd(complexAdd(TK1, complexMultiply(complex(2),TK2)), complexAdd(TK4, complexMultiply(complex(2),TK3))),complex(6)))
			// }
			// tau[s] = T;
		}

		const numerator = complexExp(complexMultiply(complexMultiply(j, tau[s]), complex(chi)));
		const denominator = complexSubtract(complexMultiply(complex(2), tau[s]), qsqr);
		// console.log(numerator, denominator);
		zs = complexAdd(zs, complexDivide(numerator, denominator));
		if (complexAbs(complexSubtract(test, zs)) < residue_options.precision * complexAbs(test)) {
			// console.log("Precision Break",test, zs);
			break;
		}
		test = zs;
	}
	// console.log(tau);

	return complexAbs(zs) * Math.fround(Math.sqrt(2*pi*chi));
}


function taustep(TA, dell, delde) {
	var a = complexPow(dell,2);
	a = complexMultiply(complex(2), complexMultiply(TA, a));
	a = complexSubtract(a, complex(1));
	a = complexDivide(delde, a);

	var b = complexPow(dell,2);
	b = complexMultiply(complex(2), complexMultiply(TA, b));
	b = complexSubtract(b, complex(1));
	b = complexDivide(complexPow(delde,2), complexPow(b,3));

	var c = complexMultiply(dell, complexPow(TA,2));
	c = complexMultiply(complex(4), c);
	c = complexAdd(c,complex(1));
	c = complexMultiply(dell, c);
	c = complexSubtract(complexMultiply(complex(2), TA), c);
	c = complexMultiply(c, dell);

	var output = complexAdd(TA, a);
	output = complexAdd(output, complexMultiply(b,c));
	return output;
}

function complexPolar(r,a) {
	return [
		r * Math.cos(a),
		r * Math.sin(a)
	]
}

function complex(a) {
	return [a,0];
}

function complexAbs([a,b]) {
	return Math.sqrt(a*a + b*b);
}

function complexSqrt([c,d]) {
	const a = Math.sqrt( (c + Math.sqrt(c*c + d*d)) / 2);
	const b = ((d<0)?-1:1) * Math.sqrt( (-c + Math.sqrt(c*c + d*d)) / 2);
	return [a,b];
}

function complexExp([a,b]) {
	const A = Math.fround(Math.exp(a) * Math.cos(b));
	const B = Math.fround(Math.exp(a) * Math.sin(b));
	return [A,B];
}

function complexCos(z) {
	return [
		Math.cos(z[0]) * Math.cosh(z[1]),
		-Math.sin(z[0])*Math.sinh(z[1])
	];
}

function complexMultiply(z1, z2) {
	return [
		z1[0] * z2[0] - z1[1] * z2[1],
		z1[0] * z2[1] + z1[1] * z2[0]
	]
}

function complexDivide([a,b], [c,d]) {
	const A = (a*c + b*d) / (c*c + d*d);
	const B = (b*c - a*d) / (c*c + d*d);
	return [A,B];
}

function complexAdd(z1, z2) {
	return [
		z1[0] + z2[0],
		z1[1] + z2[1]
	];
}

function complexSubtract(z1, z2) {
	return [
		z1[0] - z2[0],
		z1[1] - z2[1]
	];
}

function complexPow([a,b],n) {
	const r = complexAbs([a,b]);
	const theta = Math.atan2(b,a);
	const R = Math.pow(r,n);
	const T = theta * n;
	return [
		R * Math.cos(T),
		R * Math.sin(T)
	]
}

function gwFit(freq, dist, rad, field, resolution) {
	var error = 0;
	if (freq > 2) {
		freq /= 1000;
	}

	var output = {
		result: "error",
		iterations: 0
	}

	if (freq > 2 || freq < 0.5) {
		return {result: "error", detail: "Invalid Frequency"};
	}

	if (field==0) return {result: "error", detail: "Invalid measured field"};

	var input = {
		frequency: freq,
		sigma: 10,
		epsilon: 15,
		distance: dist,
		erp: 1,
		field: rad
	};

	while(Math.abs(error-1) > resolution) {
		output.iterations++;
		error = field / gwave(input);
		input.sigma *= error;

		if (input.sigma <= 0.1) {
			output.sigma = 0.1;
			output.detail = "Sigma <0.1";
			return output;
		};
		if (input.sigma > 5000) {
			output.sigma = 5000;
			output.detail = "Sigma >5000";
			return output;
		};
	}

	output.result="Success";
	output.sigma = input.sigma;
	output.detail = input;
	return output;
}

function gwaveDistance(input, resolution) {
	var error = 0;
	input.distance = 200;
	resolution = resolution ?? 0.0001;
	while (Math.abs(error - 1) > resolution) {
		error = input.target / gwave(input);
		input.distance /= error;
		// console.log(input.distance, error);
		if (input.distance < 0.1) return 0.1;
		if (input.distance > 5000) return 5000;
	}
	return input.distance;
}


function gw_chart(frequency) {
	var chart = [];
	am_sigmas.forEach( (sigma) => {
		let curve = {sigma: sigma, fields: []};
		am_distances.forEach( (dist) => {
			const input = {
				sigma: sigma,
				epsilon: (sigma==5000)?80:15,
				frequency: frequency,
				erp: 1,
				field: 100,
				distance: dist
			};
			curve.fields.push({dist: dist, field: gwave(input).output});
		});
        chart.push(curve);
	});
	return chart;
}

module.exports =  {gw_chart};