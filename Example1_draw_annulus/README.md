**VTK Example 1: Draw annulus**   <br>
David Butterworth, Oct 2018

<table style="border:none;" border="0">
<tr>
<td style="vertical-align:top;" valign="top" width="65%">Draw an annulus (disk or donut shape) using 4 different methods:  <br>
 - vtkDiskSource  <br>
 - vtkContourTriangulator  <br>
 - vtkDelaunay2D  <br>
 - Lines  <br>
  <br>
 You can choose to render:  <br>
 - the outline (shown colored white)  <br>
 - the mesh edges (shown colored blue)  <br>
 - fill with solid color  <br>
 - fill with gradient based on distance along the radius  <br>
 - fill with gradient based on angle around the circumference (a color wheel)  <br>

</td><td><img width=350 src="https://raw.githubusercontent.com/dbworth/VTK_Examples/master/Example1_draw_annulus/screenshots/screenshot1-demo.png">
</td>
</tr>
</table>

Screenshots:

<table style="border:none;" border="0">
<tr>
<td>
Solid color fill:<br>
<img width=200 src="https://raw.githubusercontent.com/dbworth/VTK_Examples/master/Example1_draw_annulus/screenshots/screenshot2-solid-fill.png">
</td>
<td>
&nbsp;
</td>
<td>
&nbsp;
</td>
</tr>

<tr>
<td>
Mesh edges (low res):<br>
<img width=200 src="https://raw.githubusercontent.com/dbworth/VTK_Examples/master/Example1_draw_annulus/screenshots/screenshot3-mesh-low-resolution.png">
</td>
<td>
Mesh edges:<br>
<img width=200 src="https://raw.githubusercontent.com/dbworth/VTK_Examples/master/Example1_draw_annulus/screenshots/screenshot4-mesh-medium-resolution.png">
</td>
<td>
Mesh edges:<br>
<img width=200 src="https://raw.githubusercontent.com/dbworth/VTK_Examples/master/Example1_draw_annulus/screenshots/screenshot5-mesh-higher-resolution.png">
</td>
</tr>

<tr>
<td>
Solid fill with mesh edges:<br>
<img width=200 src="https://raw.githubusercontent.com/dbworth/VTK_Examples/master/Example1_draw_annulus/screenshots/screenshot6-solid-with-mesh.png">
</td>
<td>
&nbsp;
</td>
<td>
&nbsp;
</td>
</tr>

<tr>
<td>
Mesh with gradient fill:<br>
<img width=200 src="https://raw.githubusercontent.com/dbworth/VTK_Examples/master/Example1_draw_annulus/screenshots/screenshot7-mesh-with-gradient-fill-radius.png">
</td>
<td>
Gradient fill (radius):<br>
<img width=200 src="https://raw.githubusercontent.com/dbworth/VTK_Examples/master/Example1_draw_annulus/screenshots/screenshot8-gradient-fill-radius.png">
</td>
<td>
Gradient fill (angle):<br>
<img width=200 src="https://raw.githubusercontent.com/dbworth/VTK_Examples/master/Example1_draw_annulus/screenshots/screenshot9-gradient-fill-angle.png">
</td>
</tr>
</table>

