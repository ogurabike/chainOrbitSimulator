
var myChart;

//テーブルの更新
function F1_setVal() {
  let cp = parseFloat(document.getElementById("F1_label_cp").innerText);  //チェーンピッチ

  let crt = parseFloat(document.getElementById("F1_select_crt").value);   //チェーンリング歯数
  let d1 = crt * cp / Math.PI;
  setInnerText(d1,"F1_label_d1","0.00");                                  //チェーンリングピッチ円直径

  let cst = parseFloat(document.getElementById("F1_select_cst").value);   //スプロケット歯数
  let d2 = cst * cp / Math.PI;
  setInnerText(d2,"F1_label_d2","0.00");                                  //スプロケットピッチ円直径

  let tpt = parseFloat(document.getElementById("F1_select_tpt").value);   //テンションプーリー歯数
  let d3 = tpt * cp / Math.PI;                                            //テンションプーリー円直径
  setInnerText(d3,"F1_label_d3","0.00");

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

//軌道計算ボタンクリック
function startCalcAndOpti() {

  document.getElementById("calcResultTbl").style.display="table-row-group";

  let returnMsg;

  let obt1 = new orbit;
  F1_getVal(obt1); // 軌道計算ボタンクリックを最初に押したときはtrue、２回目以降はfalse
  if (obt1.pulleyPositionCheck(1) == false) {
    returnMsg = "	d'0中心間最小距離(チェーンステー⇔チェーンピン)が小さすぎます。上側チェーンと下側チェーンが干渉してしまいます。";
    document.getElementById("F1_label_returnMsg").innerText = returnMsg;
    return;
  }
  if (obt1.pulleyPositionCheck(2) == false) {
    returnMsg = "	d'0中心間最小距離(チェーンステー⇔チェーンピン)が大きすぎます。チェーンとプーリーは接触しません。";
    document.getElementById("F1_label_returnMsg").innerText = returnMsg;
    return;
  }

  //最小のチェーンコマ数の確認
  obt1.orbitType = 2;
  obt1.setOrbit();
  let nc_min = obt1.optimumLinkNumber();

  F1_getVal(obt1);
  obt1.orbitType = 1;
  obt1.setOrbit();

  //仮計算した軌道からチェーンのコマ数を決定
  let nc0 = obt1.x.length-1;
  let nc = obt1.optimumLinkNumber();

  // 軌道計算ボタンクリックを最初に押したときはfalse、２回目以降はtrue
  let calcFlg;
  if (document.getElementById("F1_label_returnMsg").innerText == "") {
    calcFlg = false;
  } else {
    calcFlg = true;
  }

  // 計算２回目以降かつ仮計算軌道のチェーンコマ数に修正が不要な場合は
  // 前回計算からパラメータ変更なしみなし、これ以上の再計算は行わない。
  if ((nc == nc0) && calcFlg) {
    return;
  } else {
    if (nc < nc_min) {
      nc = nc_min;
    }

    //最適化されたチェーンのコマ数をもとに、テンションプーリーの位置を最適化(ここでobt1のチェーン軌道が最終決定される)
    let d00 = obt1.d0;
    if ((obt1.set_d0Min(nc)) && !(obt1.d0 == d00)) {
      returnMsg = "\n\nなお、d'0中心間最小距離(チェーンステー⇔チェーンピン)は、以下の通り最適化されています。\n"
        + varFormat(d00 - obt1.r3,"0.0") + "mm→" + varFormat(obt1.d0 - obt1.r3,"0.0") + "mm";

      //d'0中心間最小距離(チェーンステー⇔チェーンピン)を再設定
      document.getElementById("F1_input_dd0").value = obt1.d0 - obt1.r3;  // r3 = d3/2
      F1_setVal();  // onchangeイベントを強制実行
    } else {
      returnMsg = "最適化に失敗しました。\n	d'0中心間最小距離(チェーンステー⇔チェーンピン)を見直してください。";
      document.getElementById("F1_label_returnMsg").innerText = returnMsg;
      return;
    }
    
    //吸収可能なチェーンステー増加量を計算
    //チェーンがパツパツになるまでチェーンステー長を伸ばした状態のシミュレーション
    let obt2 = new orbit;
    F1_getVal(obt2);
    obt2.setOrbit_dlcsMax(nc);

    //obt1と0bt2のチェーンステー長の差が0.1mm以下になった時は最初からパツパツだったこととする。
    let dlcs0 = obt1.dlcs;
    if (obt2.dlcs - obt1.dlcs < 0.1 ) {
      obt1 = obt2;
      obt1.dlcs = dlcs0;
    }

    returnMsg = "上記セッティングにより、最大" + varFormat(obt2.dlcs - obt1.dlcs,"0.0") 
      + "mmのチェーンステー長の増加が吸収可能です。"
      + returnMsg;

    //計算結果表示
    setResultTable(obt1,obt2);
    
    //軌道描画
    drawChart(obt1,obt2);

    //チェーンテンショナースイング角算出
    setDeflection();

    document.getElementById("F1_label_returnMsg").innerText = returnMsg;
  }
}

//計算結果表示
function setResultTable(obt1,obt2) {
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

  setInnerText(obt2.dlcs,"F1_label_dlcs_sim2","0.0");
  setInnerText(nc_2,"F1_label_nc_sim2","0");
  setInnerText(dist_2,"F1_label_dist_sim2","0.0");
  setInnerText(obt2.d0-obt2.r3,"F1_label_dd0_sim2","0.0");
  setInnerText(obt2.d0,"F1_label_d0_sim2","0.0");
  setInnerText(swangle_2,"F1_label_swangle_sim2","0.0");

  setInnerText(swangle,"F1_label_swangle_sim1_2","0.0");
  setInnerText(swangle_2,"F1_label_swangle_sim2_2","0.0");
  setInnerText(Math.abs(swangle-swangle_2),"F1_label_delta_swangle","0.0");
}

//軌道描画
function drawChart(obt1,obt2) {
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

  //document.addEventListener('touchstart', function() {}, {passive: true});
}

//チェーンテンショナースイング角算出
function setDeflection() {
  let delta_swangle = parseFloat(document.getElementById("hidden_F1_label_delta_swangle").value);      //チェーンステー長
  let rs = parseFloat(document.getElementById("F1_input_rs").value);
  setInnerText(rs * delta_swangle * Math.PI / 180,"F1_label_deflection","0.0");
}

//テンションスプリング固定位置　初期値ボタンクリック
function set_rsIni() {
  document.getElementById("F1_input_rs").value=56.0;
  setDeflection();
}

//
function setInnerText(xval,id,format) {
  document.getElementById(id).innerText = varFormat(xval,format);
  document.getElementById("hidden_" + id).value = xval;
}