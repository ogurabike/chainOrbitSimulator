
var myChart;

//軌道計算ボタンクリック
function startCalcAndOpti() {

  if (! validation()) {
    return false;
  }

  document.getElementById("calcResultTbl").style.display="table-row-group";

  let returnMsg;

  let obt1 = new orbit;
  F1_getVal(obt1); // 軌道計算ボタンクリックを最初に押したときはtrue、２回目以降はfalse
  if (obt1.pulleyPositionCheck(1) == false) {
    returnMsg = "	d'0中心間最小距離(チェーンステー⇔チェーンピン)が小さすぎます。上側チェーンと下側チェーンが干渉してしまいます。";
    document.getElementById("F1_label_returnMsg").innerText = returnMsg;
    changeBackgroundColor("F1_cell_returnMsg",true);
    return;
  }
  if (obt1.pulleyPositionCheck(2) == false) {
    returnMsg = "	d'0中心間最小距離(チェーンステー⇔チェーンピン)が大きすぎます。チェーンとプーリーは接触しません。";
    document.getElementById("F1_label_returnMsg").innerText = returnMsg;
    changeBackgroundColor("F1_cell_returnMsg",true);
    return;
  }

  //最小のチェーンコマ数の確認
  F1_getVal(obt1);
  obt1.orbitType = 2;
  obt1.setOrbit();
  let nc_min = obt1.optimumLinkNumber();

  //仮計算した軌道からチェーンのコマ数を決定
  F1_getVal(obt1);
  obt1.orbitType = 1;
  obt1.setOrbit();
  let nc0 = obt1.x.length-1;
  let nc = obt1.optimumLinkNumber();

  // 軌道計算ボタンクリックを最初に押したときはラベル("F1_label_returnMsg")は空のはずなので、
  // 同ラベルが空の時はfalse、空でない時は２回目以降としtrue
  let calcFlg;
  returnMsg = document.getElementById("F1_label_returnMsg").innerText;
  if (returnMsg == "") {
    calcFlg = false;
  } else {
    calcFlg = true;
  }

  // 計算２回目以降かつ仮計算軌道のチェーンコマ数に修正が不要な場合は前回計算からパラメータ変更なしでの再計算とみなす。
  // パラメータ変更なしでの再計算かつSim2のチェーンステー長が設定されていない場合、これ以上の再計算は行わない。
  let dlcs_sim2 = parseFloat(document.getElementById("F1_input_dlcs_sim2").value); 
  if ((nc == nc0) && calcFlg && (dlcs_sim2 == 0)) {
    changeBackgroundColor("F1_cell_deflection_sim2",false);
    return;
  } else {
    if (nc < nc_min) {
      nc = nc_min;
    }

    //吸収可能なチェーンステー増加量を計算
    //チェーンがパツパツになるまでチェーンステー長を伸ばした状態のシミュレーション
    let obt3 = new orbit;
    F1_getVal(obt3);
    obt3.setOrbit();
    obt3.setOrbit_dlcsMax(nc);

    //obt1とobt3のチェーンステー長の差が0.1mm以下になった時は最初からパツパツだったこととする。
    let dlcs0 = obt1.dlcs;
    if (obt3.dlcs - obt1.dlcs < 0.1 ) {
      obt1 = obt3;
      obt1.dlcs = dlcs0;
    }

    let obt2 = new orbit;
    if (dlcs_sim2>0) {
      F1_getVal(obt2);
      obt2.orbitType = 1;
      obt2.dlcs=dlcs_sim2;
      if (obt2.set_d0Min(nc)) {
        if (obt2.orbitType==2) {
          obt2=obt3;
        }
      }
      changeBackgroundColor("F1_cell_deflection_sim2",true);
    } else {
      changeBackgroundColor("F1_cell_deflection_sim2",false);
    }

    //前回計算からパラメータ変更なしでの再計算の場合は結果コメントは変更しない
    //そうでない場合、中心間最小距離(チェーンステー⇔チェーンピン)最適化を行い、
    //問題なく最適化され、中心間最小距離(チェーンステー⇔チェーンピン)が変わっていれば結果コメント作成。
    if ((nc == nc0) && calcFlg ) {
      changeBackgroundColor("F1_cell_returnMsg",false);
    } else {
      let d00 = obt1.d0;

      //中心間最小距離(チェーンステー⇔チェーンピン)最適化
      if ((obt1.set_d0Min(nc)) && !(obt1.d0 == d00)) {
        returnMsg = "\n\nなお、d'0中心間最小距離(チェーンステー⇔チェーンピン)は、以下の通り最適化されています。\n"
          + varFormat(d00 - obt1.r3,"0.0") + "mm→" + varFormat(obt1.d0 - obt1.r3,"0.0") + "mm";

        returnMsg = "上記セッティングにより、最大" + varFormat(obt3.dlcs - obt1.dlcs,"0.0") 
          + "mmのチェーンステー長の増加が吸収可能です。"
          + returnMsg;

        //d'0中心間最小距離(チェーンステー⇔チェーンピン)を再設定
        document.getElementById("F1_input_dd0").value = obt1.d0 - obt1.r3;  // r3 = d3/2
        F1_setVal();  // onchangeイベントを強制実行
        changeBackgroundColor("F1_cell_returnMsg",true);

      } else {
        returnMsg = "最適化に失敗しました。\n	d'0中心間最小距離(チェーンステー⇔チェーンピン)を見直してください。";
        changeBackgroundColor("F1_cell_returnMsg",true);
        
      }
    }

    //計算結果表示
    setResultTable(obt1,obt2,obt3);
    
    //軌道描画
    drawChart(obt1,obt2,obt3);

    //チェーンテンショナースイング角算出
    setDeflection();

    document.getElementById("F1_label_returnMsg").innerText = returnMsg;
  }
}

//テンションスプリング固定位置　初期値ボタンクリック
function set_rsIni() {
  document.getElementById("F1_input_rs").value=56.0;
  setDeflection();
}

//データ入力部の値が変更された時の処理
function F1_changeVal() {
  F1_setVal();
  document.getElementById("calcResultTbl").style.display="none";
}

//テーブルの更新
function F1_setVal() {
  let cp = parseFloat(document.getElementById("F1_label_cp").innerText);  //チェーンピッチ

  let crt = parseFloat(document.getElementById("F1_select_crt").value);   //チェーンリング歯数
  let d1 = crt * cp / Math.PI;
  setInnerText(d1,"F1_label_d1","0.00");                                  //チェーンリングピッチ円直径

  let cst = parseFloat(document.getElementById("F1_select_cst").value);   //スプロケット歯数
  let d2 = cst * cp / Math.PI;
  setInnerText(d2,"F1_label_d2","0.00");                                  //スプロケットピッチ円直径

  setInnerText(crt/cst,"F1_label_gearratio","0.00");                      //ギア比計算

  let tpt = parseFloat(document.getElementById("F1_select_tpt").value);   //テンションプーリー歯数
  let d3 = tpt * cp / Math.PI;                                            //テンションプーリー円直径
  setInnerText(d3,"F1_label_d3","0.00");

  let r0 = document.getElementById("F1_input_r0").value;
  if (r0=="") {
    r0=0.0;
  } else {
    r0 = parseFloat(r0);
  }
  let r0min = (d1+d3)/2 + 5;
  document.getElementById("F1_Annotation_label_r0").innerText 
    = varFormat(r0min,"0.00") + "以上を設定できます。";
  if (r0 < r0min) {
    document.getElementById("F1_input_r0").value=varFormat(r0min,"0.00");
  }

  let dd0 = parseFloat(document.getElementById("F1_input_dd0").value);    //部品中心間最小距離(チェーンステー⇔チェーンピン)
  let d0 = dd0 + (d3 / 2);
  setInnerText(d0,"F1_label_d0","0.00");                                  //部品中心間最小距離(チェーンステー⇔プーリー)
}

//テーブル値の取得
function F1_getVal(obt) {
  obt.cp = parseFloat(document.getElementById("F1_label_cp").innerText);          //チェーンピッチ
  obt.r0 = parseFloat(document.getElementById("F1_input_r0").value);              //チェーンテンショナースイング半径
  obt.r1 = parseFloat(document.getElementById("hidden_F1_label_d1").value) / 2;   //チェーンリングピッチ円半径
  obt.r2 = parseFloat(document.getElementById("hidden_F1_label_d2").value) / 2;   //スプロケットピッチ円半径
  obt.r3 = parseFloat(document.getElementById("hidden_F1_label_d3").value) / 2;   //テンションプーリーピッチ円半径
  obt.d0 = parseFloat(document.getElementById("hidden_F1_label_d0").value);		    //中心間最小距離(チェーンステー⇔プーリー)
  obt.dlcs = parseFloat(document.getElementById("F1_input_dlcs").value);          //チェーンステー長

  //半コマチェーン
  let elements = document.getElementsByName("F1_input_chainTypes");
  let len = elements.length;
  let checkValue = "";

  for (let i = 0; i < len; i++){
      if (elements.item(i).checked){
          checkValue = elements.item(i).value;
      }
  }

  if (checkValue == "標準") {
    obt.isHalfChain = false;
  } else if (checkValue == "半コマ") {
    obt.isHalfChain = true;
  }
}

//計算結果表示
function setResultTable(obt1,obt2,obt3) {
  let nc = obt1.x.length - 1;
  let dist = Math.sqrt((obt1.x[1] - obt1.x[0]) ** 2 + (obt1.y[1] - obt1.y[0]) ** 2);
  let swangle = vbDegrees(Math.atan(obt1.d0 / Math.sqrt(obt1.r0 ** 2 - obt1.d0 ** 2)));

  setInnerText(obt1.dlcs,"F1_label_dlcs_sim1","0.0");
  setInnerText(nc,"F1_label_nc_sim1","0");
  setInnerText(dist,"F1_label_dist_sim1","0.0");
  setInnerText(obt1.d0-obt1.r3,"F1_label_dd0_sim1","0.0");
  setInnerText(obt1.d0,"F1_label_d0_sim1","0.0");
  setInnerText(swangle,"F1_label_swangle_sim1","0.0");
  
  let nc_2 = obt2.x.length - 1;
  let dist_2 = Math.sqrt((obt2.x[1] - obt2.x[0]) ** 2 + (obt2.y[1] - obt2.y[0]) ** 2);
  let swangle_2 = vbDegrees(Math.atan(obt2.d0 / Math.sqrt(obt2.r0 ** 2 - obt2.d0 ** 2)));
  
  if (obt2.dlcs=="") {
    setValue(obt2.dlcs,"F1_input_dlcs_sim2","0.0");
  }
  setInnerText(nc_2,"F1_label_nc_sim2","0");
  setInnerText(dist_2,"F1_label_dist_sim2","0.0");
  setInnerText(obt2.d0-obt2.r3,"F1_label_dd0_sim2","0.0");
  setInnerText(obt2.d0,"F1_label_d0_sim2","0.0");
  setInnerText(swangle_2,"F1_label_swangle_sim2","0.0");

  document.getElementById("F1_Annotation_label_dlcs_sim2").innerText = 
    "(" + varFormat(obt1.dlcs,"0.0")
    + "～"
    + varFormat(obt3.dlcs,"0.0") + "の間で設定できます)";

  let nc_3 = obt3.x.length - 1;
  let dist_3 = Math.sqrt((obt3.x[1] - obt3.x[0]) ** 2 + (obt3.y[1] - obt3.y[0]) ** 2);
  let swangle_3 = vbDegrees(Math.atan(obt3.d0 / Math.sqrt(obt3.r0 ** 2 - obt3.d0 ** 2)));

  setInnerText(obt3.dlcs,"F1_label_dlcs_sim3","0.0");
  setInnerText(nc_3,"F1_label_nc_sim3","0");
  setInnerText(dist_3,"F1_label_dist_sim3","0.0");
  setInnerText(obt3.d0-obt3.r3,"F1_label_dd0_sim3","0.0");
  setInnerText(obt3.d0,"F1_label_d0_sim3","0.0");
  setInnerText(swangle_3,"F1_label_swangle_sim3","0.0");

  setInnerText(Math.abs(swangle-swangle_2),"F1_label_delta_swangle_sim2","0.0");
  setInnerText(Math.abs(swangle-swangle_3),"F1_label_delta_swangle_sim3","0.0");
}

//軌道描画
function drawChart(obt1,obt2,obt3) {
  //すでにグラフが存在すれば消す
  if (myChart){
    myChart.destroy();
  }

  //チャート作成
  let ret1 = [];
  for (let mmm = 0; mmm < obt1.x.length; mmm++) {
    ret1[mmm] = { x: obt1.x[mmm], y: obt1.y[mmm] };
  }

  let ret2 = [];
  for (let mmm = 0; mmm < obt2.x.length; mmm++) {
    ret2[mmm] = { x: obt2.x[mmm], y: obt2.y[mmm] };
  }

  let ret3 = [];
  for (let mmm = 0; mmm < obt3.x.length; mmm++) {
    ret3[mmm] = { x: obt3.x[mmm], y: obt3.y[mmm] };
  }

  var ctx = document.getElementById('mychart-scatter');
  myChart = new Chart(ctx, {
    type: 'bubble',
    data: {
      datasets: [{
        label: 'Sim-1',
        data: ret1,
        backgroundColor: '#f88',
      },
      {
        label: 'Sim-2',
        data: ret2,
        backgroundColor: '#96f',
      },
      {
        label: 'Sim-3',
        data: ret3,
        backgroundColor: '#f00',
      }],
    },
    options: {
      scales: {
        yAxes: [{ ticks: {min: -200, max: 200 ,stepSize: 100}}],
        xAxes: [{ ticks: {min: -700, max: 200 ,stepSize: 100}}]
      }
    },
  });
}

//チェーンテンショナースイング角算出
function setDeflection() {
  let rs = parseFloat(document.getElementById("F1_input_rs").value);

  let delta_swangle2 = parseFloat(document.getElementById("hidden_F1_label_delta_swangle_sim2").value);      //チェーンステー長
  setInnerText(rs * delta_swangle2 * Math.PI / 180,"F1_label_deflection_sim2","0.0");

  let delta_swangle3 = parseFloat(document.getElementById("hidden_F1_label_delta_swangle_sim3").value);      //チェーンステー長
  setInnerText(rs * delta_swangle3 * Math.PI / 180,"F1_label_deflection_sim3","0.0");
}

//
function setInnerText(xval,id,format) {
  document.getElementById(id).innerText = varFormat(xval,format);
  document.getElementById("hidden_" + id).value = xval;
}

//
function setValue(xval,id,format) {
  document.getElementById(id).value = varFormat(xval,format);
  document.getElementById("hidden_" + id).value = xval;
}

function changeBackgroundColor(id,flg) {
  if (flg) {
    document.getElementById(id).style.backgroundColor = "#eef374";
  } else {
    document.getElementById(id).style.backgroundColor = "";
  }
}

//入力チェック
function validation() {
  
  if (numCheck("F1_input_dlcs","データ入力/バイク仕様/フレーム/チェーンステー長",true)
    && numCheck("F1_input_r0","データ入力/チェーンテンショナー設定/テンションプーリー取付位置/スイング半径",true)
    && numCheck("F1_input_dd0","データ入力/チェーンテンショナー設定/テンションプーリー取付位置/部品間距離 チェーン⇔チェーンステー",true)
    && numCheck("F1_input_rs","テンションスプリング仕様/テンションスプリング固定位置",true)
    && numCheck("F1_input_dlcs_sim2","テンションスプリング仕様/Sim-2/チェーンステー長",false)
    ) {
    
    let val = document.getElementById("F1_input_dlcs_sim2").value;                                //テンションスプリング仕様/Sim-2/チェーンステー長
    let lowerlimtVal = parseFloat(document.getElementById("hidden_F1_label_dlcs_sim1").value);    //テンションスプリング仕様/Sim-2/チェーンステー長
    let upperlimtVal = parseFloat(document.getElementById("hidden_F1_label_dlcs_sim3").value);    //テンションスプリング仕様/Sim-2/チェーンステー長
    if (!(val == "")) {
      val=parseFloat(val);
      if (!upperlimitCheck("テンションスプリング仕様/Sim-2/チェーンステー長",upperlimtVal,val)
        || !lowerlimitCheck("テンションスプリング仕様/Sim-2/チェーンステー長",lowerlimtVal,val)) {
          return false;
      }
    }
    return true;
  } else {
    return false;
  }
}

//必須＆半角数字入力チェック
function numCheck(id,name,isRequired ) {
  let val = document.getElementById(id).value;
  let regex = new RegExp(/\d+[\.]?\d*/);

  if (val == "") {
    if (isRequired) {
      alert('ERROR : 必須項目に値がセットされていません。\n'+ name);
      return false;
    } else {
      return true;
    }
  }

  if (regex.test(val)) {
    return true;
  } else {
    alert('ERROR : 半角数字を入力してください。\n'+ name + '\n入力値 : '+val);
    return false;
  }
}

//上限チェック
function upperlimitCheck(name,limitVal,testVal) {
  if (testVal>limitVal) {
    alert('ERROR : 入力値が上限値を超えました。\n'+ name + '\n上限値 : '+limitVal + '\n入力値 : '+testVal);
    return false;
  } else {
    return true;
  }
}

//下限チェック
function lowerlimitCheck(name,limitVal,testVal) {
  if (testVal<limitVal) {
    alert('ERROR : 入力値が下限値を下回りました。\n'+ name + '\n下限値 : '+limitVal + '\n入力値 : '+testVal);
    return false;
  } else {
    return true;
  }
}

//指定桁数で丸め込み
function varFormat(val1, str1) {
  let val2 = str1.length - 2;
  if (val2 == -1) { val2 = 0 }
  return Math.round(val1 * 10 ** val2) / (10 ** val2);
}

//rad→degree変換
function vbDegrees(val1) {
  return 180 / Math.PI * val1;
}