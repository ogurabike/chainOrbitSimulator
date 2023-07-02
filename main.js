class orbit {
  isHalfChain;
  cp;
  r0;
  r1;
  r2;
  r3;
  d0;
  dlcs;
  px;
  py;
  x;
  y;

  constructor() {
    this.cp = parseFloat(document.getElementById("F1_label_cp").innerText);      //チェーンピッチ
    this.r0 = parseFloat(document.getElementById("F1_input_r0").value);          //チェーンテンショナースイング半径
    this.r1 = parseFloat(document.getElementById("F1_label_d1").innerText) / 2;  //チェーンリングピッチ円半径
    this.r2 = parseFloat(document.getElementById("F1_label_d2").innerText) / 2;  //スプロケットピッチ円半径
    this.r3 = parseFloat(document.getElementById("F1_label_d3").innerText) / 2;  //テンションプーリーピッチ円半径
    this.d0 = parseFloat(document.getElementById("F1_label_d0").innerText);			  //中心間最小距離(チェーンステー⇔プーリー)
    this.dlcs = parseFloat(document.getElementById("F1_input_dlcs").value);      //チェーンステー長

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
      this.isHalfChain = false;
    } else if (checkValue == "半コマ") {
      this.isHalfChain = true;
    }

    this.px = [];
    this.py = [];
    this.x = [];
    this.y = [];

    for (let mmm = 0; mmm <= 6; mmm++) {
      this.px[mmm] = 0;
      this.py[mmm] = 0;
    }
  }

  //チェーン軌道の計算_chainOrbit
  setOrbit() {
    this.x = [];
    this.y = [];
    this.px = [];
    this.py = [];

    let anglex;
    let angle1;
    let angle2;

    // P6,P1の座標計算
    this.px[1] = -this.r1 * (this.r1 - this.r2) / this.dlcs;
    this.py[1] = this.r1 * Math.sqrt(this.dlcs ** 2 - (this.r1 - this.r2) ** 2) / this.dlcs;

    this.px[6] = -(this.dlcs ** 2 + this.r2 * (this.r1 - this.r2)) / this.dlcs;
    this.py[6] = this.r2 * Math.sqrt(this.dlcs ** 2 - (this.r1 - this.r2) ** 2) / this.dlcs;

    // P2,P3の座標計算
    let delta;
    let deltamin;

    delta = 0;
    deltamin = 0;
    angle1 = 0;
    anglex = Math.PI + Math.atan(this.d0 / Math.sqrt(this.r0 ** 2 - this.d0 ** 2));

    for (let mmm = 1; mmm <= 100; mmm++) {
      if (this.r1 == this.r2) {
        anglex = anglex
          + mmm / 100 * (Math.PI / 2 - Math.atan(this.d0 / Math.sqrt(this.r0 ** 2 - this.d0 ** 2)));
      } else {
        anglex = anglex
          + mmm / 100 * (
            Math.atan(Math.sqrt(this.dlcs ** 2 - (this.r1 - this.r2) ** 2) / (this.r1 - this.r2))
            - Math.atan(this.d0 / Math.sqrt(this.r0 ** 2 - this.d0 ** 2))
          );
      }

      delta = Math.abs(this.r3 - Math.abs(Math.sqrt(this.r0 ** 2 - this.d0 ** 2) * Math.cos(anglex) + this.d0 * Math.sin(anglex) + this.r1));

      if (mmm == 1) {
        deltamin = delta;
        angle1 = anglex;
      } else if (delta < deltamin) {
        deltamin = delta;
        angle1 = anglex;
      } else if (mmm > 1 && delta > deltamin) {
        break;
      }

    }

    this.px[2] = this.r1 * Math.cos(angle1);
    this.py[2] = this.r1 * Math.sin(angle1);

    this.px[3] = -Math.sqrt(this.r0 ** 2 - this.d0 ** 2) + this.r3 * Math.abs(Math.cos(angle1));
    this.py[3] = -this.d0 + this.r3 * Math.abs(Math.sin(angle1));

    // P4,P5の座標計算
    delta = 0;
    deltamin = 0;
    angle2 = 0;

    for (let mmm = 1; mmm <= 100; mmm++) {
      anglex = (1 + (1 / 100) * mmm) * Math.PI;
      delta = Math.abs(
        this.r3 - Math.abs((this.dlcs - Math.sqrt(this.r0 ** 2 - this.d0 ** 2)) * Math.cos(anglex) - this.d0 * Math.sin(anglex) - this.r2)
      );

      if (mmm == 1) {
        deltamin = delta;
        angle2 = anglex;
      } else if (delta < deltamin) {
        deltamin = delta;
        angle2 = anglex;
      } else if (mmm > 1
        && Math.ceil(delta * 10) > Math.ceil(deltamin * 10)) {

        //P4のY座標はC3のY座標より大きい
        if (-this.d0 - this.r3 * Math.sin(angle2) < -this.d0) {
          deltamin = Math.abs(this.r3 - Math.abs((this.dlcs - Math.sqrt(this.r0 ** 2 - this.d0 ** 2)) * Math.cos(Math.PI) - this.d0 * Math.sin(Math.PI) - this.r2));
        } else {
          break;
        }
      }
    }

    this.px[5] = this.r2 * Math.cos(angle2) - this.dlcs;
    this.py[5] = this.r2 * Math.sin(angle2);

    this.px[4] = -Math.sqrt(this.r0 ** 2 - this.d0 ** 2) + this.r3 * Math.cos(angle2);
    this.py[4] = -this.d0 - this.r3 * Math.sin(angle2);

    //P1→P2

    this.x[1] = this.px[1];
    this.y[1] = this.py[1];

    setOnCircumference(0, 0, this.r1, this.px[1], this.py[1], this.px[2], this.py[2], this.x, this.y, this.cp, true);

    //P2→P3
    tangent(this.px[2], this.py[2], this.px[3], this.py[3], this.x, this.y, this.cp);

    //P3→P4
    setOnCircumference(
      -Math.sqrt(this.r0 ** 2 - this.d0 ** 2), -this.d0, this.r3, this.px[3], this.py[3], this.px[4], this.py[4], this.x, this.y, this.cp, false);

    //P4→P5
    tangent(this.px[4], this.py[4], this.px[5], this.py[5], this.x, this.y, this.cp);

    //P5→P6
    setOnCircumference(-this.dlcs, 0, this.r2, this.px[5], this.py[5], this.px[6], this.py[6], this.x, this.y, this.cp, true);

    //P6→P1
    tangent(this.px[6], this.py[6], this.px[1], this.py[1], this.x, this.y, this.cp);

    this.x[0] = this.x[this.x.length - 1];
    this.y[0] = this.y[this.y.length - 1];
  }

  //最長チェーンステー長計算_Optimisation2
  setOrbit_dlcsMax() {

    //　チェーンコマ数取得
    let nc;
    this.setOrbit();

    let dcp0 = 0;
    let dcp1 = 0;

    nc = this.optimumLinkNumber();

    // チェーンステー長導出
    let dlt = 1.0;

    for (let nnn = 0; nnn <= 10000; nnn++) {
      dcp0 = Math.abs((this.x[1] - this.x[0]) ** 2 + (this.y[1] - this.y[0]) ** 2 - this.cp ** 2);
      this.dlcs = this.dlcs + dlt;

      this.px = [];
      this.py = [];
      this.x = [];
      this.y = [];

      for (let mmm = 0; mmm <= 6; mmm++) {
        this.px[mmm] = 0;
        this.py[mmm] = 0;
      }

      this.set_dlcsMax();

      dcp1 = Math.abs((this.x[1] - this.x[0]) ** 2 + (this.y[1] - this.y[0]) ** 2 - this.cp ** 2);

      // console.log("nnn:"+nnn+" nc:"+nc+" this.x.length - 1:"+(this.x.length - 1)+" dcp0:"+dcp0+" dcp1:"+dcp1+" dlcs:"+this.dlcs+" dlt:"+dlt);

      if ((dlt == 1.0) && (nc+1 == this.x.length - 1)) {
        this.dlcs = this.dlcs - 4*dlt;
        dlt = 0.1;
      } else if ((dlt == 0.1) && (nc+1 == this.x.length - 1)) {
        this.dlcs = this.dlcs - 4*dlt;
        dlt = 0.05;
      } else if ((dlt == 0.05) && (nc+1 == this.x.length - 1)) {
        this.dlcs = this.dlcs - 4*dlt;
        dlt = 0.01;
      } else if ((dlt == 0.01) && (nc+1 == this.x.length - 1)) {
        this.dlcs = this.dlcs - 4*dlt;
        dlt = 0.001;
      } else if ((dlt == 0.001) && (nc == this.x.length - 1)
        && (Math.abs((this.x[1] - this.x[0]) ** 2 + (this.y[1] - this.y[0]) ** 2 - this.cp ** 2) < 0.1)) {
        break;
      } else if (nc < this.x.length - 1) {
        return false;
      }
    }
    return true;
  }
  
  // チェーンコマ数の最適化
  optimumLinkNumber() {
    let nc = this.x.length - 1;

    // 最適化不要時の処理
    if (this.isHalfChain) {
      if (Math.abs((this.x[1] - this.x[0]) ** 2 + (this.y[1] - this.y[0]) ** 2 - this.cp ** 2) <= 0.1) {
        return nc;
      }
    } else if (Math.ceil(nc / 2) == nc / 2
      && Math.abs((this.x[1] - this.x[0]) ** 2 + (this.y[1] - this.y[0]) ** 2 - this.cp ** 2) <= 0.1) {
      return nc;
    }

    // コマ数調整
    if (this.isHalfChain) {
      nc = nc - 1;   //半コマチェーンの時は１コマ詰める
    } else {
      if (Math.ceil(nc / 2) == nc / 2) {
        nc = nc - 2; //標準ピッチでコマ数が偶数の時は2コマ詰める
      } else {
        nc = nc - 1; //標準ピッチでコマ数が偶数の時は1コマ詰める
      }
    }
    return nc;
  }
  
  // 中心間最小距離(チェーンステー⇔チェーンピン)最適化_Optimisation
  set_d0Min() {

    let nc;
    let nc2;

    this.setOrbit();

    nc = this.x.length - 1;
    nc2 = this.optimumLinkNumber();

    if (nc==nc2) {
      return true;
    }

    let dlt = 10.0;
    for (let mmm = 1; mmm <= 100000; mmm++) {
      this.d0 = this.d0 + dlt;
      this.px = [];
      this.py = [];
      this.x = [];
      this.y = [];
      nc = 0;

      if (this.isContact() == false) {
        return false;
      }

      this.setOrbit();
      nc = this.x.length - 1;

      if ((dlt == 10.0) && (nc <= nc2)) {
        this.d0 = this.d0 - dlt;
        dlt = 1.0;
      } else if ((dlt == 1.0) && (nc <= nc2)) {
        this.d0 = this.d0 - dlt;
        dlt = 0.1;
      } else if ((dlt == 0.1) && (nc <= nc2)) {
        this.d0 = this.d0 - dlt;
        dlt = 0.01;
      } else if ((dlt == 0.01) && (nc <= nc2)) {
        this.d0 = this.d0 - dlt;
        dlt = 0.001;
      } else if ((dlt == 0.001) && (Math.abs((this.x[1] - this.x[0]) ** 2 + (this.y[1] - this.y[0]) ** 2 - this.cp ** 2) < 0.1)) {
        break;
      } else if ((dlt == 0.001) && (nc < nc2)) {
        return 0;	//最適化失敗
      }
    }
    return true;
  }

  //チェーン軌道の計算（テンションプーリ無し）_chainOrbit2
  set_dlcsMax() {
    this.x = [];
    this.y = [];
    this.px = [];
    this.py = [];

    // P6,P1の座標計算
    this.px[1] = -this.r1 * (this.r1 - this.r2) / this.dlcs;
    this.py[1] = this.r1 * Math.sqrt(this.dlcs ** 2 - (this.r1 - this.r2) ** 2) / this.dlcs;

    this.px[6] = -(this.dlcs ** 2 + this.r2 * (this.r1 - this.r2)) / this.dlcs;
    this.py[6] = this.r2 * Math.sqrt(this.dlcs ** 2 - (this.r1 - this.r2) ** 2) / this.dlcs;

    // P2,P5の座標計算
    this.px[2] = this.px[1];
    this.py[2] = -this.py[1];

    this.px[5] = this.px[6];
    this.py[5] = -this.py[6];

    // P1→P2
    this.x[1] = this.px[1];
    this.y[1] = this.py[1];

    setOnCircumference(0, 0, this.r1, this.px[1], this.py[1], this.px[2], this.py[2], this.x, this.y, this.cp, true);

    // P2→P5
    tangent(this.px[2], this.py[2], this.px[5], this.py[5], this.x, this.y, this.cp);

    // P5→P6
    setOnCircumference(-this.dlcs, 0, this.r2, this.px[5], this.py[5], this.px[6], this.py[6], this.x, this.y, this.cp, true);

    // P6→P1
    tangent(this.px[6], this.py[6], this.px[1], this.py[1], this.x, this.y, this.cp);

    this.x[0] = this.x[this.x.length - 1];
    this.y[0] = this.y[this.y.length - 1];

    // d0導出
    let dlt = 10.0;
    for (let nnn = 0; nnn <= 20000; nnn++) {
      this.d0 = this.d0 + dlt;
      if ((this.isContact() == false) && (dlt == 10.0)) {
        this.d0 = this.d0 - dlt;
        dlt = 1.0;
      } else if ((this.isContact() == false) && (dlt == 1.0)) {
        this.d0 = this.d0 - dlt;
        dlt = 0.1;
      } else if ((this.isContact() == false) && (dlt == 0.1)) {
        this.d0 = this.d0 - dlt;
        dlt = 0.01;
      } else if ((this.isContact() == false) && (dlt == 0.01)) {
        this.d0 = this.d0 - dlt;
        break;
      }
    }
    return true;
  }
  
  //プーリー⇔チェーン接触チェック
  //d0を増やしすぎて、接線P6P1をチェーンステーの軸中心(y=0)で反転させた直線より
  //下側にプーリーが行ってしまうとチェーンとプーリーが接しなくなってしまう。
  isContact() {

    let slope;
    let intercept;
    let p1x;
    let p1y;
    let p6x;
    let p6y;
    let c3x;
    let c3y;

    p1x = -this.r1 * (this.r1 - this.r2) / this.dlcs;
    p1y = this.r1 * Math.sqrt(this.dlcs ** 2 - (this.r1 - this.r2) ** 2) / this.dlcs;

    p6x = -(this.dlcs ** 2 + this.r2 * (this.r1 - this.r2)) / this.dlcs;
    p6y = this.r2 * Math.sqrt(this.dlcs ** 2 - (this.r1 - this.r2) ** 2) / this.dlcs;

    slope = linearfunction(p1x, -p1y, p6x, -p6y,)["slope"];
    intercept = linearfunction(p1x, -p1y, p6x, -p6y,)["intercept"];

    if (this.r0 < this.d0) {
      return false;
    }

    c3x = -Math.sqrt(this.r0 ** 2 - this.d0 ** 2);
    c3y = -this.d0;

    if ((slope * c3x - c3y + intercept) / Math.sqrt(slope ** 2 + 1) <= this.r3) {
      return true;
    } else {
      return false;
    }
  }
}

function F1_setVal() {
  let cp = parseFloat(document.getElementById("F1_label_cp").innerText);  //チェーンピッチ

  let crt = parseFloat(document.getElementById("F1_select_crt").value);   //チェーンリング歯数
  document.getElementById("F1_label_d1").innerText = crt * cp / Math.PI;        //チェーンリングピッチ円直径

  let cst = parseFloat(document.getElementById("F1_select_cst").value);   //スプロケット歯数
  document.getElementById("F1_label_d2").innerText = cst * cp / Math.PI;        //スプロケットピッチ円直径

  let tpt = parseFloat(document.getElementById("F1_select_tpt").value);   //テンションプーリー歯数
  let d3 = tpt * cp / Math.PI;                                                //テンションプーリー円直径
  document.getElementById("F1_label_d3").innerText = d3;

  let dd0 = parseFloat(document.getElementById("F1_input_dd0").value);    //部品中心間最小距離(チェーンステー⇔チェーンピン)
  document.getElementById("F1_label_d0").innerText = dd0 + (d3 / 2);            //部品中心間最小距離(チェーンステー⇔プーリー)
}

//いきなり最適化を実施
function startCalcAndOpti() {
  let dd0 = parseFloat(document.getElementById("F1_input_dd0").value);    //部品中心間最小距離(チェーンステー⇔チェーンピン)
  //startCalc();
  startOpti();
  document.getElementById("F1_input_dd0").value =dd0;
}

//計算開始ボタンクリック
function startCalc() {

  F1_setVal();

  let mmm;
  let nnn;

  let returnMsg;

  let obt1 = new orbit;

  if (obt1.isContact() == false) {
    returnMsg = "	d'0中心間最小距離(チェーンステー⇔チェーンピン)が大きすぎます。チェーンとプーリーは接触しません。";
    return returnMsg;
  }

  obt1.setOrbit();
  
  let nc = obt1.x.length - 1;
  let dist = vbFormat(Math.sqrt((obt1.x[1] - obt1.x[0]) ** 2 + (obt1.y[1] - obt1.y[0]) ** 2), "0.00");
  let swangle = vbFormat(vbDegrees(Math.atan(obt1.d0 / Math.sqrt(obt1.r0 ** 2 - obt1.d0 ** 2))), "0.000");

  document.getElementById("F1_label_dlcs").innerText = obt1.dlcs;
  document.getElementById("F1_label_nc").innerText = nc;
  document.getElementById("F1_label_dist").innerText = dist;
  document.getElementById("F1_label_dd0_1").innerText = obt1.d0-obt1.r3;
  document.getElementById("F1_label_d0_1").innerText = obt1.d0;
  document.getElementById("F1_label_swangle").innerText = swangle;

  let obt2 = new orbit;
  obt2.setOrbit_dlcsMax();
  
  let nc_2 = obt2.x.length - 1;
  let dist_2 = vbFormat(Math.sqrt((obt2.x[1] - obt2.x[0]) ** 2 + (obt2.y[1] - obt2.y[0]) ** 2), "0.00");
  let swangle_2 = vbFormat(vbDegrees(Math.atan(obt2.d0 / Math.sqrt(obt2.r0 ** 2 - obt2.d0 ** 2))), "0.000");

  document.getElementById("F1_label_maxdlcs").innerText = obt2.dlcs;
  document.getElementById("F1_label_nc_2").innerText = nc_2;
  document.getElementById("F1_label_dist_2").innerText = dist_2;
  document.getElementById("F1_label_dd0_2").innerText = obt2.d0-obt2.r3;
  document.getElementById("F1_label_d0_2").innerText = obt2.d0;
  document.getElementById("F1_label_swangle_2").innerText = swangle_2;

  document.getElementById("F1_label_swangle__2").innerText = swangle;
  document.getElementById("F1_label_swangle_2__2").innerText = swangle_2;
  document.getElementById("F1_label_delta_swangle").innerText = Math.abs(swangle-swangle_2);

  setDeflection();

  let ret1 = [];
  for (let mmm = 0; mmm < obt1.x.length; mmm++) {
    ret1[mmm] = { x: obt1.x[mmm], y: obt1.y[mmm] };
  }

  let ret2 = [];
  for (let mmm = 0; mmm < obt2.x.length; mmm++) {
    ret2[mmm] = { x: obt2.x[mmm], y: obt2.y[mmm] };
  }

  var ctx = document.getElementById('mychart-scatter');
  var myChart = new Chart(ctx, {
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
        y: { min: -150, max: 150 },
        x: { min: -700, max: 200 },
      },
    },
  });

  returnMsg = "チェーン軌道の計算が完了しました。\n"
    + "\n";

  //半コマチェーン
  if (obt1.isHalfChain) {
    if (Math.abs((obt1.x[1] - obt1.x[0]) ** 2 + (obt1.y[1] - obt1.y[0]) ** 2 - obt1.cp ** 2) <= 0.1) {
      returnMsg = returnMsg + "最適化は不要です。\n開始コマと終端コマの距離がチェーンピッチ(" + obt1.cp + "mm)と一致しています。";
    } else {
      returnMsg = returnMsg + "最適化が必要です。\n開始コマと終端コマの距離がチェーンピッチ(" + obt1.cp + "mm)と一致していないため、チェーンの先端と終端が接続できません。";
    }
  }

  //標準ピッチ
  else {

    //標準ピッチかつコマ数が偶数
    if (Math.ceil((obt1.x.length - 1) / 2) == (obt1.x.length - 1) / 2) {
      if (Math.abs((obt1.x[1] - obt1.x[0]) ** 2 + (obt1.y[1] - obt1.y[0]) ** 2 - obt1.cp ** 2) <= 0.1) {
        returnMsg = returnMsg + "最適化は不要です。\n開始コマと終端コマの距離がチェーンピッチ(" + obt1.cp + "mm)と一致しています。";
      } else {
        returnMsg = returnMsg + "最適化が必要です。\n開始コマと終端コマの距離がチェーンピッチ(" + obt1.cp + "mm)と一致していないため、\nチェーンの先端と終端が接続できません。";
      }
    }

    //標準ピッチなのにコマ数が奇数
    else {
      returnMsg = returnMsg + "最適化が必要です。\nコマ数が奇数になっているため、チェーンの先端と終端が接続できません\n（チェーンタイプ：標準選択時は偶数のみ）。";
    }
  }

  document.getElementById("F1_label_returnMsg").innerText = returnMsg;

  return returnMsg;
}

//最適化ボタンクリック
function startOpti() {
  let returnMsg;
  let obt1 = new orbit;
  obt1.setOrbit();

  let d00 = obt1.d0;

  obt1.set_d0Min();

  if (obt1.d0 == d00) {

    returnMsg = "最適化の必要はありません。\n"
      + "チェーンコマ数=" + (obt1.x.length - 1) + "pcs.\n"
      + "チェーン接合部（開始コマと終端コマ）コマ間距離="
      + vbFormat(Math.sqrt((obt1.x[1] - obt1.x[0]) ** 2 + (obt1.y[1] - obt1.y[0]) ** 2), "0.00") + " mm";

  } else if (obt1.d0 == 0) {

    returnMsg = "最適化に失敗しました。\n	d'0中心間最小距離(チェーンステー⇔チェーンピン)を見直してください。";

  } else {

    returnMsg = "最適化が完了しました。\n d'0中心間最小距離(チェーンステー⇔チェーンピン)が変更されています。\n" + (d00 - obt1.r3) + "mm→" + (obt1.d0 - obt1.r3) + "mm";
    document.getElementById("F1_input_dd0").value = obt1.d0 - obt1.r3; // r3 = d3/2
    startCalc();

  }
  document.getElementById("F1_label_returnMsg").innerText = returnMsg;
  return returnMsg;
}

//円周上に指定ピッチで点を配置するプログラム
function setOnCircumference(
  x0, y0,
  r,
  x1, y1,
  x2, y2,
  x, y,
  d,
  isClockwise) {
  let s;
  let t;
  let u;
  let v;

  let rad;
  let angle;
  let dd;
  let ddd;

  let rad1;
  let rad2;
  let rad12;
  let deltarad12;
  let ndeltarad12;

  let nc;

  //円弧（円周）の両端の点から配置できる点の最大数を求める。
  if (x1 >= x0) {
    if (y1 >= y0) {
      //第1象限
      rad1 = Math.acos(Math.abs(x1 - x0) / r);
    } else {
      //第4象限
      rad1 = 2 * Math.PI - Math.acos(Math.abs(x1 - x0) / r);
    }
  } else {
    if (y1 >= y0) {
      //第2象限
      rad1 = Math.PI - Math.acos(Math.abs(x1 - x0) / r);
    } else {
      //第3象限
      rad1 = Math.PI + Math.acos(Math.abs(x1 - x0) / r);
    }
  }

  if (x2 >= x0) {
    if (y2 >= y0) {
      //第1象限
      rad2 = Math.acos(Math.abs(x2 - x0) / r);
    } else {
      //第4象限
      rad2 = 2 * Math.PI - Math.acos(Math.abs(x2 - x0) / r);
    }
  } else {
    if (y2 >= y0) {
      //第2象限
      rad2 = Math.PI - Math.acos(Math.abs(x2 - x0) / r);
    } else {
      //第3象限
      rad2 = Math.PI + Math.acos(Math.abs(x2 - x0) / r);
    }
  }

  if (isClockwise) {
    if (rad1 < rad2) {
      rad12 = 2 * Math.PI - (rad2 - rad1);
    } else {
      rad12 = rad1 - rad2;
    }
  } else {
    if (rad1 < rad2) {
      rad12 = rad2 - rad1;
    } else {
      rad12 = 2 * Math.PI - (rad1 - rad2);
    }
  }

  do {
    if (rad12 > 2 * Math.PI) {
      rad12 = rad12 - 2 * Math.PI;
    }
  } while (rad12 > 2 * Math.PI);

  deltarad12 = Math.abs(2 * Math.asin(d / 2 / r));
  ndeltarad12 = Math.floor(rad12 / deltarad12);

  //接線上の点(s,t)から指定距離dの位置にある円周上の最初の点(u,v)を求める
  nc = x.length - 1;
  s = x[nc];
  t = y[nc];
  ddd = 100000;
  for (let mmm = 1; mmm <= 20000; mmm++) {
    if (isClockwise) {
      rad = rad1 - (mmm / 20000) * deltarad12;
    } else {
      rad = rad1 + (mmm / 20000) * deltarad12;
    }

    u = r * Math.cos(rad) + x0;
    v = r * Math.sin(rad) + y0;
    dd = Math.abs((s - u) ** 2 + (t - v) ** 2 - d ** 2);

    if (ddd > dd) {
      ddd = dd;
      angle = rad;
    }
  }

  nc = nc + 1
  x[nc] = r * Math.cos(angle) + x0;
  y[nc] = r * Math.sin(angle) + y0;

  //残りの点を配置する。
  for (let mmm = 1; mmm <= ndeltarad12 + 5; mmm++) {

    if (isClockwise) {
      angle = angle - deltarad12;
      if (rad1 > rad2 && angle < rad2) {
        break;
      } else if (rad1 < rad2 && angle < rad2 - 2 * Math.PI) {
        break;
      }
    } else {
      angle = angle + deltarad12;
      if (rad1 < rad2 && angle > rad2) {
        break;
      } else if (rad1 > rad2 && angle > rad2 + 2 * Math.PI) {
        break;
      }
    }

    nc = nc + 1;
    x[nc] = r * Math.cos(angle) + x0;
    y[nc] = r * Math.sin(angle) + y0;
  }

}

//両端を指定された接線上に指定ピッチで点を配置するプログラム
function tangent(
  p1x, p1y,
  p2x, p2y,
  x, y,
  d) {

  let nc;
  let s;
  let t;
  let u;
  let v;

  let a;
  let b;
  let c;

  let slope;
  let intercept;

  let angle;

  //両端の点から接線の式を求める
  slope = linearfunction(p1x, p1y, p2x, p2y)["slope"];
  intercept = linearfunction(p1x, p1y, p2x, p2y)["intercept"];

  //円周上の点(s,t)から指定距離dの位置にある接線上の最初の点(u,v)を求める
  nc = x.length - 1;
  s = x[nc];
  t = y[nc];

  a = 1 + slope ** 2;
  b = 2 * (slope * (intercept - t) - s);
  c = s ** 2 + (intercept - t) ** 2 - d ** 2;

  u = solv2(a, b, c)["x1"];
  v = slope * u + intercept;

  if (isOnLineSegment(p1x, p1y, p2x, p2y, u, v) == false) {

    u = solv2(a, b, c)["x2"];
    v = slope * u + intercept;

    if (isOnLineSegment(p1x, p1y, p2x, p2y, u, v) == false) {
      u = 0;
      v = 0;
    }

  }

  nc = nc + 1;
  x[nc] = u;
  y[nc] = v;

  //残りの点を配置する
  angle = Math.atan(slope);

  do {
    nc = nc + 1;

    if (slope == 0) {
      x[nc] = x[nc - 1] + d;
      y[nc] = y[nc - 1];
    } else {
      x[nc] = x[nc - 1] + (p2x - p1x) / Math.abs((p2x - p1x)) * d * Math.abs(Math.cos(angle));
      y[nc] = y[nc - 1] + (p2y - p1y) / Math.abs((p2y - p1y)) * d * Math.abs(Math.sin(angle));
    }
  } while (isOnLineSegment(p1x, p1y, p2x, p2y, x[nc], y[nc]));

  //最後の点は線分からはみ出してるから削除する
  x.pop();
  y.pop();
}

//二次関数における解と係数の関係
function solv2(a, b, c) {
  let ret = [];
  let g;
  let h;

  g = -b / 2 / a
  h = Math.sqrt(b ** 2 - 4 * a * c) / 2 / a

  ret["x1"] = g - h
  ret["x2"] = g + h
  return ret;
}

//二点を通る直線の傾きと切片の計算
function linearfunction(
  x1, y1,
  x2, y2
) {

  let ret = [];
  ret["slope"] = (y1 - y2) / (x1 - x2);
  ret["intercept"] = -ret["slope"] * x1 + y1;

  return ret;

}

//線分上の点であることチェック
function isOnLineSegment(
  x1, y1,
  x2, y2,
  x, y) {

  let dist0;
  let dist1;
  let dist2;

  //線分(x1,y1)→(x2,y2)の長さ
  dist0 = Math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2);

  //点(x1,y1)→中点(x,y)の距離
  dist1 = Math.sqrt((x - x1) ** 2 + (y - y1) ** 2);

  //点(x2,y2)→中点(x,y)の距離
  dist2 = Math.sqrt((x - x2) ** 2 + (y - y2) ** 2);

  if (vbFormat(dist1 + dist2, "0.000") == vbFormat(dist0, "0.000")) {
    return true;
  } else {
    return false;
  }

}

//指定桁数で丸め込み
function vbFormat(val1, str1) {
  let val2 = str1.length - 2;

  if (val2 == -1) { val2 = 0 }
  return Math.round(val1 * 10 ** val2) / (10 ** val2);
}

//rad→degree変換
function vbDegrees(val1) {
  return 180 / Math.PI * val1;
}

function setDeflection() {
  let delta_swangle = parseFloat(document.getElementById("F1_label_delta_swangle").innerText);      //チェーンステー長
  let rs = parseFloat(document.getElementById("F1_input_rs").value);
  document.getElementById("F1_label_deflection").innerText = rs * delta_swangle * Math.PI / 180;
}

function set_rsIni() {
  document.getElementById("F1_input_rs").value=56.0;
  setDeflection();
}