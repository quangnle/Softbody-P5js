let W = 800, H = 600;
let nPoints = 50;
let Radius = Math.min(190.0, W / 2 - 1);
let R2 = Radius * Radius;
let BallRadius = 0.516;

//Spring constants
let KS             = 755.0;    //755 by default
let KD             = 35.0;     //35.0 by default
let FAPP           = 110.0;

let Gx = 0;
let Gy = 110.0;
let GMagnitude = 110.0;
let Mass = 1;
let FinalPressure = 30000;

let deltaTime = 0.01;

let TangentialDampingFactor = 0.99;                         
let NormalDampingFactor = 0.1;

let points = [];
let springs = [];

let pressure = 0;
let mouseP = false;

let ball = new SoftBall(0.516, 50, 30000, Radius);

function setup(){
	var canvas = createCanvas(W, H);
	canvas.parent('sketch-holder');
	ball.create();
}

function draw() {
	background(255);
	stroke(0);
	fill(255);
	ellipse(W/2, H/2, 2*Radius, 2*Radius);
	
	ball.draw(mouseX, mouseY);
}

function mousePressed() {
	ball.mouseClicked = true;
}

function mouseReleased() {
	ball.mouseClicked = false;
}

