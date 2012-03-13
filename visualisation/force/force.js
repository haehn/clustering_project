var width = 1024,
    height = 800,
    node,
    link,
    root;
    

var color = d3.scale.category20();


var svg = d3.select("#chart").append("svg:svg")
    .attr("width", width)
    .attr("height", height)
    .attr("pointer-events", "all")
    .call(d3.behavior.zoom().on("zoom", redraw))
    .on("mouseup",function() {
        shiftingPoint = false}
     )
    .append('svg:g');

svg.append('svg:rect')
    .attr('width', width)
    .attr('height', height)
    .attr('fill', 'rgba(1,1,1,0)');


var shiftingPoint = false;
function redraw() { 
    // If zooming, change transformation 
    if (shiftingPoint || d3.event.sourceEvent.target.nodeName == "circle") {
        shiftingPoint = true;
        return;
    }
    else {
        console.log("here", d3.event.translate, d3.event.scale, d3.event.sourceEvent.target.nodeName );
        svg.attr("transform", "translate(" + d3.event.translate + ")scale(" + d3.event.scale + ")");
    } 
}

var force = d3.layout.force()
    .gravity(0.01)
    .charge(function(d){ return( d.isLeaf ? -500 : -100) })
    .friction(.5)
    .linkDistance(function(d) { return 10*d.length; })
    .linkStrength(function(d) { 
         return d.width; 
     })
    .size([width, height]);

d3.json("graph.json", function(json) {
  root = json;
  link = svg.selectAll("line.link")
      .data(json.links)
      .enter().append("line")
      .attr("class", "link")
      .style("stroke-width", function(d) { return d.width*10; })
      .on("click",function(d,i) {
           selectClasses = intersectSorted(d.source.class,d.target.class);
           colorClasses(selectClasses);
       });

  var legend = d3.select('#legend').selectAll("span.legendItem")
      .data(json.legend)
      .enter().append("span")
      .text(function(d,i) { return("Class "+(i+1))})
      .attr("class", "button")
      .attr("classid",function(d,i){return(i+1)})
      .on("click",function(d,i) { 
           colorClasses([i+1]);
       });

  node = svg.selectAll("g.node")
      .data(json.nodes)
      .enter().append("svg:g")
      .attr("class", "node")
      .call(force.drag)
      .on("mouseup",function() {
        shiftingPoint = false}
       );


   node.append("svg:text")
        .attr("class", "nodetext")
        .attr("dx", 10)
        .attr("dy", ".35em")
        .text(function(d) { return(d.isLeaf ? d.name : "")});

   node.append("svg:circle")
      .attr("class","circle")
      .attr("r", 5)
      .style("fill", function(d) { return(d.isLeaf ? color(1) : 'lightgrey') });

   //node.append("title")
   //   .text(function(d) { return d.name; });



  force
      .nodes(json.nodes)
      .links(json.links)
      .on("tick", function(e) {

         link.attr("x1", function(d) { return d.source.x; })
             .attr("y1", function(d) { return d.source.y; })
             .attr("x2", function(d) { return d.target.x; })
             .attr("y2", function(d) { return d.target.y; });
     
         node.attr("cx", function(d) { return d.x; })
             .attr("cy", function(d) { return d.y; });
        // node.attr("transform", function(d) { return "translate(" + d.x + "," + d.y + ")"; });

       })
     .start();

   function colorClasses(k) {
       d3.selectAll("span.buttondown").attr("class", "button");
       d3.selectAll("span.button").filter(function(d2,i2) {
           return(intersectSorted([i2+1],k).length > 0 ? true : false);
       }).attr("class", "buttondown");

       d3.selectAll("line.link").style("stroke","grey").attr("selected",false);
       d3.selectAll("line.link").filter(function(d2) {
           tmp = intersectSorted(d2.source.class,d2.target.class);
           return(intersectSorted(tmp,k).length > 0 ? true : false);
       }).style("stroke","red").attr("selected",true);

   }
});

// the following function comes from:
//    http://stackoverflow.com/a/1885660
function intersectSorted(a, b)
{
  var ai=0, bi=0;
  var result = new Array();

  while( ai < a.length && bi < b.length )
  {
     if      (a[ai] < b[bi] ){ ai++; }
     else if (a[ai] > b[bi] ){ bi++; }
     else /* they're equal */
     {
       result.push(a[ai]);
       ai++;
       bi++;
     }
  }

  return result;
}
function minusSorted(a, b)
{
  var ai=0, bi=0;
  var result = new Array();

  while( ai < a.length && bi < b.length )
  {
     if      (a[ai] < b[bi] ){ 
        result.push(a[ai]);
        ai++; 
     }
     else if (a[ai] > b[bi] ){ bi++; }
     else /* they're equal */
     {
       ai++;
       bi++;
     }
  }

  return result;
}
