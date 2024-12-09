<template>
    <div class='background'>
        <div class='center-div'>
            <hr>
            <el-row id='header-col1' > <!-- Title of page -->
                <el-col  :span="2">
                    <span class='hid'>"      "</span>
                </el-col>
                <el-col  :span="19">
                       <span class='hbig'>4D V</span>
                       <span class='hsmall'>IRTUAL </span>
                       <span class='hbig'>R</span>
                       <span class='hsmall'>EGENERATING </span>
                       <span class='hbig'>H</span>
                       <span class='hsmall'>EART</span>
                </el-col>
                <el-col :span="2">
                    <img src="./assets/HF2.png" class="logo1" >
                </el-col>
            </el-row>
            <el-tabs v-model="activeName" type="border-card" > <!-- The main menu -->
                <el-tab-pane label="Atlas" name="atlas"> <!-- The Atlas panel -->
                    <el-row> <!-- section header -->
                        <el-col :span="19"> 
                            <h2>The Spatiotemporal Atlas of Regenerating Zebrafish Heart</h2>
                        </el-col>
                        <el-col :span="5"> 
                            <h2>The 3D atlas of zebrafish heart</h2>
                        </el-col>
                    </el-row>
                    <el-row> <!-- The second row of database -->
                        <el-col :span="19" > 
                            <div class='box'>
                            <el-row> <!-- The second row of database -->
                                <el-col :span="6">  <!-- The E14-16h card-->
                                    <div class='card'>
                                        <el-row> <!-- card name -->
                                            <el-col :span="24">
                                                <h3>uninjured</h3>
                                            </el-col>
                                        </el-row>
                                        <el-row> <!-- card image -->
                                            <el-col :span="24">
                                                <img src="./assets/T0.png" class="image" >
                                            </el-col>
                                        </el-row>
                                        <el-row> <!-- card buttons -->
                                            <el-col :span="24">
                                                <h3><a href="/zebrafish_VRH/heart3d_ctrl_3D/index.html?atlas=/zebrafish_VRH/heart3d_ctrl_3D/" target="_blank">Click to view in 3D</a></h3>
                                            </el-col>
                                        </el-row>
                                    </div>
                                </el-col>
                                <el-col :span="6">  <!-- The E14-16h card-->
                                    <div class='card'>
                                        <el-row> <!-- card name -->
                                            <el-col :span="24">
                                                <h3>6 hpa</h3>
                                            </el-col>
                                        </el-row>
                                        <el-row> <!-- card image -->
                                            <el-col :span="24">
                                                <img src="./assets/T2.png" class="image" >
                                            </el-col>
                                        </el-row>
                                        <el-row> <!-- card buttons -->
                                            <el-col :span="24">
                                                <h3><a href="/zebrafish_VRH/heart3d_6h_3D/index.html?atlas=/zebrafish_VRH/heart3d_6h_3D/" target="_blank">Click to view in 3D</a></h3>
                                            </el-col>
                                        </el-row>
                                    </div>
                                </el-col>
                                <el-col :span="6">  <!-- The E14-16h card-->
                                    <div class='card'>
                                        <el-row> <!-- card name -->
                                            <el-col :span="24">
                                                <h3>12 hpa</h3>
                                            </el-col>
                                        </el-row>
                                        <el-row> <!-- card image -->
                                            <el-col :span="24">
                                                <img src="./assets/T3.png" class="image" >
                                            </el-col>
                                        </el-row>
                                        <el-row> <!-- card buttons -->
                                            <el-col :span="24">
                                                <h3><a href="/zebrafish_VRH/heart3d_12h_3D/index.html?atlas=/zebrafish_VRH/heart3d_12h_3D/" target="_blank">Click to view in 3D</a></h3>
                                            </el-col>
                                        </el-row>
                                    </div>
                                </el-col>
                                <el-col :span="6">  <!-- The E14-16h card-->
                                    <div class='card'>
                                        <el-row> <!-- card name -->
                                            <el-col :span="24">
                                                <h3>1 dpa</h3>
                                            </el-col>
                                        </el-row>
                                        <el-row> <!-- card image -->
                                            <el-col :span="24">
                                                <img src="./assets/T4.png" class="image" >
                                            </el-col>
                                        </el-row>
                                        <el-row> <!-- card buttons -->
                                            <el-col :span="24">
                                                <h3><a href="/zebrafish_VRH/heart3d_1d_3D/index.html?atlas=/zebrafish_VRH/heart3d_1d_3D/" target="_blank">Click to view in 3D</a></h3>
                                            </el-col>
                                        </el-row>
                                    </div>
                                </el-col>
                            </el-row>
                            <el-row> <!-- The third row of database -->
                                <el-col :span="6">  <!-- The E14-16h card-->
                                    <div class='card'>
                                        <el-row> <!-- card name -->
                                            <el-col :span="24">
                                                <h3>3 dpa</h3>
                                            </el-col>
                                        </el-row>
                                        <el-row> <!-- card image -->
                                            <el-col :span="24">
                                                <img src="./assets/T5.png" class="image" >
                                            </el-col>
                                        </el-row>
                                        <el-row> 
                                            <el-col :span="24">
                                                <h3><a href="/zebrafish_VRH/heart3d_3d_3D/index.html?atlas=/zebrafish_VRH/heart3d_3d_3D/" target="_blank">Click to view in 3D</a></h3>
                                            </el-col>
                                        </el-row>
                                    </div>
                                </el-col>
                                <el-col :span="6">  <!-- The E14-16h card-->
                                    <div class='card'>
                                        <el-row> <!-- card name -->
                                            <el-col :span="24">
                                                <h3>7 dpa</h3>
                                            </el-col>
                                        </el-row>
                                        <el-row> <!-- card image -->
                                            <el-col :span="24">
                                                <img src="./assets/T6.png" class="image" >
                                            </el-col>
                                        </el-row>
                                        <el-row> <!-- card buttons -->
                                            <el-col :span="24">
                                                <h3><a href="/zebrafish_VRH/heart3d_7d_3D/index.html?atlas=/zebrafish_VRH/heart3d_7d_3D/" target="_blank">Click to view in 3D</a></h3>
                                            </el-col>
                                        </el-row>
                                    </div>
                                </el-col>
                                <el-col :span="6">  
                                    <div class='card'>
                                        <el-row> <!-- card name -->
                                            <el-col :span="24">
                                                <h3>14 dpa</h3>
                                            </el-col>
                                        </el-row>
                                        <el-row> <!-- card image -->
                                            <el-col :span="24">
                                                <img src="./assets/T7.png" class="image" >
                                            </el-col>
                                        </el-row>
                                        <el-row> <!-- card buttons -->
                                            <el-col :span="24">
                                                <h3><a href="/zebrafish_VRH/heart3d_14d_3D/index.html?atlas=/zebrafish_VRH/heart3d_14d_3D/" target="_blank">Click to view in 3D</a></h3>
                                            </el-col>
                                        </el-row>
                                    </div>
                                </el-col>
                                <el-col :span="6">  <!-- The E14-16h card-->
                                    <div class='card'>
                                        <el-row> <!-- card name -->
                                            <el-col :span="24">
                                                <h3>28 dpa</h3>
                                            </el-col>
                                        </el-row>
                                        <el-row> <!-- card image -->
                                            <el-col :span="24">
                                                <img src="./assets/T8.png" class="image" >
                                            </el-col>
                                        </el-row>
                                        <el-row> <!-- card buttons -->
                                            <el-col :span="24">
                                                <h3><a href="/zebrafish_VRH/heart3d_28d_3D/index.html?atlas=/zebrafish_VRH/heart3d_28d_3D/" target="_blank">Click to view in 3D</a></h3>
                                            </el-col>
                                        </el-row>
                                    </div>
                                </el-col>
                            </el-row>
                            <el-row>
                                <h2>Note: data in this box are regeneration datasets and their 3D coordinate were obtained by registrate the 2D regeneration slices into 3D uninjured atlas.</h2>
                            </el-row>
                            </div>
                        </el-col>
                        <el-col :span="5"> 
                            <div class="box">
                            <el-row> <!-- The row of database -->
                                <el-col :span="34">  <!-- The E14-16h card-->
                                    <div class='card'>
                                        <el-row> <!-- card name -->
                                            <el-col :span="24">
                                                <h3>High-res uninjured heart</h3> <h3 style="color:red">(recommend)</h3>
                                            </el-col>
                                        </el-row>
                                        <el-row> <!-- card image -->
                                            <el-col :span="24">
                                                <img src="./assets/ST3D_merge.png" class="image" >
                                            </el-col>
                                        </el-row>
                                        <el-row> <!-- card buttons -->
                                            <el-col :span="24">
                                                <h3><a href="/zebrafish_VRH/heart3d_merge/index.html?atlas=/zebrafish_VRH/heart3d_merge/" target="_blank">Click to view</a></h3>
                                            </el-col>
                                        </el-row>
                                    </div>
                                </el-col>
                            </el-row>
                            <el-row> <!-- The row of database -->
                                <el-col :span="24"> 
                                    <div class='card'>
                                        <el-row> 
                                            <el-col :span="24">
                                                <h3>Stereo-seq raw atlas</h3>  <h3 style="color:white">"                 " </h3>
                                            </el-col>
                                        </el-row>
                                        <el-row> 
                                            <el-col :span="24">
                                                <img src="./assets/ST3D_raw.png" class="image" >
                                            </el-col>
                                        </el-row>
                                        <el-row> 
                                            <el-col :span="24">
                                                <h3><a href="/zebrafish_VRH/heart3d/index.html?atlas=/zebrafish_VRH/heart3d/" target="_blank">Click to view</a></h3>
                                            </el-col>
                                        </el-row>
                                    </div>
                                </el-col> 
                            </el-row>
                            </div>
                        </el-col>
                    </el-row>
                </el-tab-pane>
            </el-tabs>
            <hr>
            <div class='footer'>
                <el-row>
                    <el-col :span="4">
                        <h1>          </h1>
                    </el-col>
                    <el-col :span="16">
                        <img src="./assets/logos.png" class="logo" >
                    </el-col>
                    <el-col :span="4">
                        <h1>       </h1>
                    </el-col>
                </el-row>
                <el-row>
                    <el-col :span="24">
                        <p>Copyright &#169; All Rights Reserved 2023 BGI Research <a href=https://beian.miit.gov.cn/>鲁ICP备2022034375号</a></p>
                    </el-col>
                </el-row>
            </div>
        </div>
    </div>
</template>

<script>
export default {
    data() {
        return {
          tid: '1',
          activeName: 'atlas',
        }
    },
};
</script>

<style>
  #header-col1 {
    background-color: #93c572;
    text-align: center;
    line-height: 150px;
    font-weight: bold;
    color: white;
  }
  .hid{
    color: #93c572;
  }
  .hbig{
    text-align: right;
    font-size: 50px;
    font-weight: bold;
    color: white;
  }
  .hsmall{
    text-align: left;
    font-size: 40px;
    font-weight: bold;
    color: white;
  }
  .logo1{
     width:130px;
     margin:0px;
     display:block;
     padding:10px;
  } 
  .logo{
    width: 1000px;
    display: block;
    margin:10px;
    padding:20px;
  }
  .footer{
    border-color: #eeeeee;
    border-style: solid;
    border-width: 1px;
    color: #8b8c8f;
    font-size: 14px;
    padding:10px;
  }
  .image {
    width: 100%;
    display: block;
  }
  .box{
     border-width: 2px;
     border-style: solid;
     border-color: #93c572;
     margin:2px
  }
  .card{
     border-width: 2px;
     border-style: solid;
     border-color:#111111;
     margin:10px
  }
  .tutorial{
     background-color: #ccc;
     border-width: 2px;
     border-style: dashed;
     border-color:#111111;
     text-align: left;
     margin:10px;
     padding:10px;
  }
  .left_border{
     border-left: 3px dashed gray;
  }
  .upper_border{
     border-top: 3px dashed gray;
  }
  .center-div
  {
    margin: 0 auto;
    width: 1500px;
    background-color: white;
    padding:5px;
  }
  .background
  {
      background-color: #ccc;
  }
</style>
