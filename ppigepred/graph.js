var nodes = [];
var edges = [];
var hiddenEdges = [];
var gravityConstant = 0.03;
var repulsionForce = 2000;
var degreeForceConstant = 0.25;
var springForceConstant = 0.2;
var edgeAngleRepulsionContant = 0.001;
let center;
let pg;
var neighbours = new Map();

let gravityForces;
let chargeRepulsion;
let degreeRepulsion;
let springForces;
let edgeAngleRepulsion;
let minEdgeScore = 0;
let updateEdgesFlag = false;

var minDegree = 1000000;
var maxDegree = -1;
var maxScore = 1000;

var globalPosX = 0;
var globalPosY = 0;

var maxDensity;
var minDensity;

function degree(n) {
  let degree = (neighbours.get(n) || []).length;
  minDegree = min(minDegree, degree);
  maxDegree = max(maxDegree, degree);
  return degree;
}

function percentageToHsl(percentage, hue0, hue1) {
  var hue = (percentage * (hue1 - hue0)) + hue0;
  return 'hsl(' + int(hue) + ', 100%, 50%)';
}

// Create a new node with the given position, id and name.
function newNode(pos, id, name, density, isReference) {
  return {
    id: id,
    pos: pos,
    name: name,
    density: density,
    force: createVector(0, 0),
    isReference: isReference,
    mass: () => 10*degree(id)/maxDegree + 5, // normalize the mass
    update: function() {
      let velocity = p5.Vector.div(this.force, this.mass());
      this.pos.add(velocity);
    },
    draw: function() {
      noStroke();
      let id = parseInt(this.id);
      if (hoveredNode !== null &&
          neighbours != null &&
          neighbours.has(id) &&
          neighbours.get(id).indexOf(parseInt(hoveredNode.id)) != -1) {
        stroke('red');
      }
      if (this.isReference) {
        fill('black');
      } else {
        let colorValue = (density - minDensity)/maxDensity;
        let logisticColor = 1/(1 + Math.exp(-10*colorValue)) - 0.5;
        let hsl = percentageToHsl(logisticColor, 120, 0);
        fill(color(hsl));
      }
      let diameter = 5 * this.mass();
      circle(this.pos.x, this.pos.y, diameter);
    }
  };
}


function onlyUnique(value, index, array) {
  return array.indexOf(value) === index;
}

// Places the nodes at random positions.
function initNodes(graphData) {
  for (let [id, nodeData] of Object.entries(graphData.nodes)) {
      let x = random(width);
      let y = random(height);
    let node = newNode(
        createVector(x, y),
        id,
        nodeData.name,
        nodeData.density,
        nodeData.is_reference
      );
      nodes.push(node);
  }
}


// Takes two nodes and adds each one of them as a neigbour
// of the other in the global "neighbours" map.
function addNeighbours(neighbours, n1, n2) {
  if (neighbours.has(n1)) {
    neighbours.get(n1).push(n2);
  } else {
    neighbours.set(n1, [n2])
  }
  if (neighbours.has(n2)) {
    neighbours.get(n2).push(n1);
  } else {
    neighbours.set(n2, [n1])
  }
}

function initEdges(graphData) {
  for (var i = 0; i < graphData.edges.length; i++) {
    let n1 = graphData.edges[i][0];
    let n2 = graphData.edges[i][1];
    let score = graphData.edges[i][2];

    edges.push([nodes[n1], nodes[n2], score]);
    addNeighbours(neighbours, n1, n2);
  }
}

// Applies a force of gravity to all nodes which pulls them
// to the origin.
function applyGravity() {
  nodes.forEach(node => {
    let gravity = p5.Vector.sub(node.pos, center).mult(-1).mult(gravityConstant).mult(node.mass());
    node.force = gravity;
  });
}

// Applies repulsion forces between all nodes. There are two types of such forces:
// - Charge repulsion acting between any pair of nodes. The force becomes
//   stronger when the nodes are closer.
// - Degree repulsion makes high-degree nodes repulse other high-degree nodes.
function applyRepulsionForces() {
  for (let i = 0; i < nodes.length; i++) {
    let node1 = nodes[i];
    for (let j = i + 1; j < nodes.length; j++) {
      let node2 = nodes[j];
      let dir = node1.pos.copy().sub(node2.pos);
      let dist = dir.mag();

      // apply charge repulsion
      if (chargeRepulsion) {
        let force = dir.copy().div(dist * dist).mult(repulsionForce);
        node1.force.add(force);
        node2.force.sub(force);
      }

      // apply degree rupulsion
      if (degreeRepulsion) {
        let node1Deg = degree(i);
        let node2Deg = degree(j);
        let degreeForce = degreeForceConstant*((node1Deg + 1)*(node2Deg + 1))/dist;
        node1.force.add(dir.copy().mult(degreeForce));
        node2.force.sub(dir.copy().mult(degreeForce));
      }
    }
  }
}

var h = 0;

// Applies forces on nodes in attempt to keep the length of the edges
// to their optimal length (calculated based on the number of nodes).
// If the edge is too short, the nodes incident to it would be pushed away
// from each other and if the edge is too long, the nodes would be pulled
// together.
function applySpringForces() {
  edges.forEach(edge => {
    let node1 = edge[0];
    let node2 = edge[1];
    let dir = node1.pos.copy().sub(node2.pos);
    let dist = dir.mag();
    // The optimal length depends on the score of the protein interaction.
    // We also add a term that depends on the number of edges and nodes to
    // make more dense networks better spread out.
    optimalLen = maxScore / edge[2] + Math.log(edges.length + nodes.length);
    let correction = springForceConstant*(dist - optimalLen)/Math.max(optimalLen, dist);
    node1.force.sub(dir.copy().mult(correction));
    node2.force.add(dir.copy().mult(correction));
  });
}

// Applies forces so that edges incident to a node are distributed
// evenly around the node. Basically, trying to keep a 360/deg(node)
// degree angle between them.
function applyEdgeAngleRepulsion() {
  for (let centralNode of nodes) {
    let allAdjacent = neighbours.get(centralNode.id);
    if (allAdjacent == null) {
      continue;
    }
    let adjacent = allAdjacent.filter((v, i, a) => a.indexOf(v) === i);
    let optimalAngle = 2*PI/adjacent.length;
    for (let i = 0; i < adjacent.length; i++) {
      let nodeA = nodes[adjacent[i]];
      for (let j = i + 1; j < adjacent.length; j++) {
        let nodeB = nodes[adjacent[j]];
        let centralToADir = nodeA.pos.copy().sub(centralNode.pos);
        let centralToBDir = nodeB.pos.copy().sub(centralNode.pos);
        let angle = Math.abs(centralToADir.angleBetween(centralToBDir));
        let correctionAngle = edgeAngleRepulsionContant*(optimalAngle - angle);
        let newDirCentralToA = centralToADir.copy().rotate(correctionAngle);
        let desiredPosA = newDirCentralToA.copy().add(centralNode.pos);
        let force = desiredPosA.copy().sub(nodeA.pos);

        nodeA.force.add(force);
      }
    }
  };
}

// Applies all forces for a single iteration.
function applyForces() {
  if (gravityForces) {
    applyGravity();
  }

  applyRepulsionForces();

  if (springForces) {
    applySpringForces();
  }

  if (edgeAngleRepulsion) {
    applyEdgeAngleRepulsion();
  }
}

function readData() {
    var data;
    var req = new XMLHttpRequest();
    req.open("GET", '/subgraph', false);
    req.onreadystatechange = function () {
        if(req.readyState === 4) {
            if(req.status === 200 || req.status == 0) {
                data = req.responseText;
            }
        }
    }
    req.send(null);
    return JSON.parse(data);
}

function setup() {
  gravityForces = document.getElementById("gravity").checked;
  chargeRepulsion = document.getElementById("charge").checked;
  degreeRepulsion = document.getElementById("degree").checked;
  springForces = document.getElementById("spring").checked;
  edgeAngleRepulsion = document.getElementById("angle").checked;

  let bodyWidth = document.body.clientWidth;
  let bodyHeight = document.body.clientHeight;
  createCanvas(bodyWidth, bodyHeight - 100);
  center = createVector(0, 0);
  
  graphData = readData();
  initNodes(graphData);
  initEdges(graphData);
  
  maxDensity = Math.max(...nodes.map(function(node) { return node.isReference ? 0 : node.density }));
  minDensity = Math.min(...nodes.map(function(node) { return node.isReference ? 1 : node.density }));
}

// Updates the currentScale variable so that the whole
// graph fits on the canvas.
function rescale() {
  let minX = Math.min.apply(null, nodes.map(n => n.pos.x));
  let maxX = Math.max.apply(null, nodes.map(n => n.pos.x));
  let minY = Math.min.apply(null, nodes.map(n => n.pos.y));
  let maxY = Math.max.apply(null, nodes.map(n => n.pos.y));

  let x = maxX - minX;
  let y = maxY - minY;

  return 0.8 * 1/Math.max(x/width, y/height, 1);
}

function updateEdges() {
  var newEdges = [];
  var newHiddenEdges = [];
  var newNeighbours  = new Map();

  for (let edge of edges) {
    let score = edge[2];
    if (score >= minEdgeScore) {
      newEdges.push(edge);
      addNeighbours(newNeighbours, edge[0], edge[1]);
    } else {
      newHiddenEdges.push(edge);
    }
  }
  for (let edge of hiddenEdges) {
    let score = edge[2];
    if (score >= minEdgeScore) {
      newEdges.push(edge);
      addNeighbours(newNeighbours, edge[0], edge[1]);
    } else {
      newHiddenEdges.push(edge);
    }
  }
  
  edges = newEdges;
  hiddenEdges = newHiddenEdges;
  neighbours = newNeighbours;
}

function updateForces() {
  gravityConstant = 0.1;
  repulsionForce = 10*maxDegree*Math.log(nodes.length + edges.length);
  degreeForceConstant = (Math.log(maxDegree) * Math.log(nodes.length + edges.length))/5000;
  springForceConstant = Math.log(nodes.length + edges.length)/100;
  edgeAngleRepulsionContant = 0.001;
}

var currentScale = 1;
var hoveredNode = null;

function draw() {
  // white background
  background(258);

  translate(width/2, height/2);

  currentScale = rescale();
  if (currentScale != 1) {
    scale(currentScale);
  }
  
  updateForces();
 
  if (updateEdgesFlag) {
    updateEdges();
    updateEdgesFlag = false;
  }

  hoveredNode = null;
  let translatedMouseX = mouseX - width/2;
  let translatedMouseY = mouseY - height/2;
  for (let node of nodes) {
    let r = currentScale * 2.5 * node.mass();
    if (dist(node.pos.x*currentScale, node.pos.y*currentScale, translatedMouseX, translatedMouseY) < r) {
      hoveredNode = node;
    }
  }

  edges.forEach(edge => {
    stroke(0);
    if (hoveredNode !== null && (edge[0].id == hoveredNode.id || edge[1].id == hoveredNode.id)) {
      stroke('red');
    }
    line(edge[0].pos.x, edge[0].pos.y, edge[1].pos.x, edge[1].pos.y);
  });

  nodes.forEach(node => {
    node.draw();
    node.update();
  });

  if (hoveredNode !== null) {
    textSize(14/currentScale);
    fill('black');
    stroke('black');
    text(hoveredNode.name + ' (' + hoveredNode.density + ')', hoveredNode.pos.x + 10, hoveredNode.pos.y - 10)
  }

  applyForces();
}