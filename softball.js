function distance(p1, p2){
	return Math.sqrt((p2.x - p1.x)*(p2.x - p1.x) + (p2.y - p1.y)*(p2.y - p1.y));
}

var SoftBall = function(radius, nPoints, finalPressure, boundary){
	let KS = 755.0;    
	let KD = 35.0;     
	let FAPP  = 110.0;
	let Gx = 0;
	let Gy = 110.0;
	let GMagnitude = 110.0;
	let Mass = 1;
	let dt = 0.01;
	let TanDF = 0.99;                         
	let NormDF = 0.1;
	
	this.radius = radius;
	this.nPoints = nPoints;
	this.finalPressure = finalPressure;	
	this.boundary = boundary;
	this.points = [];
	this.springs = [];
	this.mouseClicked = false;
	this.pressure = 0;
	
	this.draw = function(mx, my){
		this.mx = mx;
		this.my = my;
		
		fill('#0F0');
		noStroke();
		beginShape();
		for (let i=0; i < this.points.length; i++) {
			vertex(this.points[i].x, this.points[i].y);
		}
		endShape(CLOSE);

		stroke('#F00');
		if(this.mouseClicked) {
			line(this.mx, this.my, this.points[0].x, this.points[0].y);
		}
		
		let nSteps = 2;
		let pStep = this.finalPressure / 300;
		for (let i = 0; i < nSteps; i++){
			this.accumulateForces();
			this.integrateHeun();
			if (this.pressure < this.finalPressure) {
				this.pressure += pStep;
			}
		}
	}

	this.create = function(){
		for (let i = 0; i < this.nPoints; i++) {
			let x = this.radius * Math.cos(i * 2 * Math.PI / this.nPoints) + (W >> 1);
			let y = this.radius * Math.sin(i * 2 * Math.PI / this.nPoints) + (H >> 1);
			this.points.push(new Point(x, y, 0, 0, 0, 0));
		}
		
		for (let i = 0; i < this.nPoints; i++) {
			let spring = new Spring(this.points[i], this.points[(i + 1) % this.nPoints], 0, 0, distance(this.points[i], this.points[(i + 1) % this.nPoints]));
			this.springs.push(spring);
		}
	}
	
	this.accumulateForces = function(){
		let x1, x2, y1, y2;
		let d;
		let vx12;
		let vy12;
		let f;
		let fx0, fy0;
		let volume = 0;
		let pressurev = 0;
		
		for (let i = 0; i < this.nPoints; i++) {
			this.points[i].fx = (this.pressure >= this.finalPressure) ? Gx * Mass : 0;
			this.points[i].fy = (this.pressure >= this.finalPressure) ? Gy * Mass : 0;
		}
		
		if(this.mouseClicked) {
			x1 = this.points[0].x;
			y1 = this.points[0].y;
			x2 = this.mx;
			y2 = this.my;

			d = distance(this.points[0], {x:x2, y:y2});
			f = (d - 2.2) * 10 + (this.points[0].vx * (x1 - x2) + this.points[0].vy * (y1 - y2)) / d;

			fx0 = ((x1 - x2) / d ) * f;
			fy0 = ((y1 - y2) / d ) * f;

			this.points[0].fx -= fx0;
			this.points[0].fy -= fy0;
		}

		for (let i = 0; i < this.nPoints; i++) {
			x1 = this.springs[i].p1.x;
			y1 = this.springs[i].p1.y;
			x2 = this.springs[i].p2.x;
			y2 = this.springs[i].p2.y;
			
			d = distance(this.springs[i].p1, this.springs[i].p2);
			if (d != 0) {
				vx12 = this.springs[i].p1.vx - this.springs[i].p2.vx;
				vy12 = this.springs[i].p1.vy - this.springs[i].p2.vy;

				f = (d - this.springs[i].length) * KS + (vx12 * (x1 - x2) + vy12 * (y1 - y2)) * KD / d;

				fx0 = ((x1 - x2) / d ) * f;
				fy0 = ((y1 - y2) / d ) * f;

				this.springs[i].p1.fx -= fx0;
				this.springs[i].p1.fy -= fy0;
						   
				this.springs[i].p2.fx += fx0;
				this.springs[i].p2.fy += fy0;
			}
			
			this.springs[i].nx = (y2 - y1)/d;
			this.springs[i].ny = (x1 - x2)/d;
		}
	  
		for (let i = 0; i < this.nPoints; i++) {
			x1 = this.springs[i].p1.x;
			y1 = this.springs[i].p1.y;
			x2 = this.springs[i].p2.x;
			y2 = this.springs[i].p2.y;

			d = distance(this.springs[i].p1, this.springs[i].p2);
			volume += 0.5 * Math.abs(x1 - x2) * Math.abs(this.springs[i].nx) * d;
		}
	  
		for (let i = 0; i < this.nPoints; i++) {
			x1 = this.springs[i].p1.x;
			y1 = this.springs[i].p1.y;
			x2 = this.springs[i].p2.x;
			y2 = this.springs[i].p2.y;

			d = distance(this.springs[i].p1, this.springs[i].p2);
			pressurev = d * this.pressure * (1.0/volume);

			this.springs[i].p1.fx += this.springs[i].nx * pressurev;
			this.springs[i].p1.fy += this.springs[i].ny * pressurev;
			this.springs[i].p2.fx += this.springs[i].nx * pressurev;
			this.springs[i].p2.fy += this.springs[i].ny * pressurev;
		}
	}
	
	this.integrateHeun = function(){		
		let R2 = this.boundary * this.boundary;
		let drx, dry;

		let fxsaved = [];
		let fysaved = [];
		let vxsaved = [];
		let vysaved = [];

		for (let i = 0; i < this.nPoints; i++) {
			fxsaved.push(this.points[i].fx);
			fysaved.push(this.points[i].fy);    
			vxsaved.push(this.points[i].vx);
			vysaved.push(this.points[i].vy);

			this.points[i].vx += this.points[i].fx / Mass * dt;
			drx = this.points[i].vx * dt;
			this.points[i].x += drx;

			this.points[i].vy += this.points[i].fy / Mass * dt;
			dry = this.points[i].vy * dt;
			this.points[i].y += dry;
		}

		this.accumulateForces();

		for (let i=0; i < this.nPoints; i++) {
			this.points[i].vx = vxsaved[i] + (this.points[i].fx + fxsaved[i])/Mass * dt/2;
			drx = this.points[i].vx * dt;
			this.points[i].x += drx;

			this.points[i].vy = vysaved[i] + (this.points[i].fy + fysaved[i])/Mass * dt/2;
			dry = this.points[i].vy * dt;
			this.points[i].y += dry;

			this.points[i].x = Math.min(this.points[i].x, (W>>1) + this.boundary);
			this.points[i].x = Math.max(this.points[i].x, (W>>1) - this.boundary);
			this.points[i].y = Math.min(this.points[i].y, (H>>1) + this.boundary);
			this.points[i].y = Math.max(this.points[i].y, (H>>1) - this.boundary);

			if (this.points[i].x + drx >  Math.sqrt(R2 - Math.pow(this.points[i].y - (H>>1), 2)) + (W>>1) || 
				this.points[i].x + drx < -Math.sqrt(R2 - Math.pow(this.points[i].y - (H>>1), 2)) + (W>>1)){
				drx *= -1;                           
				dry *= -1;                           

				let vx0 = this.points[i].vx;
				let vy0 = this.points[i].vy;

				let sinTheta = (this.points[i].y - H*0.5) / this.boundary;
				let cosTheta = (this.points[i].x - W*0.5) / this.boundary;
				let sincos = sinTheta * cosTheta;
				let sinsin = sinTheta * sinTheta;
				let coscos = cosTheta * cosTheta;
				

				this.points[i].vx = -vx0;
				this.points[i].vy = -vy0;
				this.points[i].vx = vy0 * (-TanDF * sincos - NormDF * sincos) + vx0 * (TanDF * sinsin - NormDF * coscos);
				this.points[i].vy = vy0 * (TanDF * coscos - NormDF * sinsin) + vx0 * (-TanDF * sincos - NormDF * sincos);
			}

			if (this.points[i].y > ((H + Radius) >> 1)) { 
				this.points[i].y = Math.min(this.points[i].y,  Math.sqrt(Math.abs(R2 - Math.pow(this.points[i].x - W*0.5, 2))) + H*0.5);
			}                                                                                                       
																													
			if (this.points[i].y < ((H - Radius) >> 1)) {                                                           
				this.points[i].y = Math.max(this.points[i].y, -Math.sqrt(Math.abs(R2 - Math.pow(this.points[i].x - W*0.5, 2))) + H*0.5);
			}                                                                                                      
																													
			if (this.points[i].x > ((W + Radius) >> 1)) {                                                          
				this.points[i].x = Math.min(this.points[i].x,  Math.sqrt(Math.abs(R2 - Math.pow(this.points[i].y - H*0.5, 2))) + W*0.5);
			}                                                                                                       
																													
			if (this.points[i].x < ((W - Radius) >> 1)) {                                                           
				this.points[i].x = Math.max(this.points[i].x, -Math.sqrt(Math.abs(R2 - Math.pow(this.points[i].y - H*0.5, 2))) + W*0.5);
			}
		}
	}
}