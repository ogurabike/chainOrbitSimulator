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
    orbitType;  //1:プーリー有り　2:プーリー無し

    constructor() {
        this.px = [];
        this.py = [];
        this.x = [];
        this.y = [];
      
        for (let mmm = 0; mmm <= 6; mmm++) {
            this.px[mmm] = 0;
            this.py[mmm] = 0;
        }
    }

    setOrbit() {
        let delta;
        let deltamin;
        let anglex;
        let angle1;
        let angle2;
        
        this.x = [];
        this.y = [];
        this.px = [];
        this.py = [];

        // P6の座標計算
        this.px[6] = -(this.dlcs ** 2 + this.r2 * (this.r1 - this.r2)) / this.dlcs;
        this.py[6] = this.r2 * Math.sqrt(this.dlcs ** 2 - (this.r1 - this.r2) ** 2) / this.dlcs;
        
        // P1の座標計算
        this.px[1] = -this.r1 * (this.r1 - this.r2) / this.dlcs;
        this.py[1] = this.r1 * Math.sqrt(this.dlcs ** 2 - (this.r1 - this.r2) ** 2) / this.dlcs;

        // P2,P3の座標計算
        if (this.orbitType == 1) {
            delta = 0;
            deltamin = 0;
            angle1 = 0;
            anglex = Math.PI + Math.abs(Math.atan(this.d0 / Math.sqrt(this.r0 ** 2 - this.d0 ** 2)));

            for (let mmm = 1; mmm <= 100; mmm++) {
                if (this.r1 == this.r2) {
                    anglex = anglex
                    + mmm / 100 * Math.abs(Math.PI / 2 - Math.atan(this.d0 / Math.sqrt(this.r0 ** 2 - this.d0 ** 2)));
                } else {
                    anglex = anglex
                    + mmm / 100 * (
                        Math.abs(Math.atan(Math.sqrt(this.dlcs ** 2 - (this.r1 - this.r2) ** 2) / (this.r1 - this.r2))
                        - Math.atan(this.d0 / Math.sqrt(this.r0 ** 2 - this.d0 ** 2)))
                    );
                    //console.log(Math.atan(Math.sqrt(this.dlcs ** 2 - (this.r1 - this.r2) ** 2) / (this.r1 - this.r2)));
                    //console.log(Math.atan(this.d0 / Math.sqrt(this.r0 ** 2 - this.d0 ** 2)));
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

        } else if (this.orbitType == 2) {
            this.px[2] = this.px[1];
            this.py[2] = -this.py[1];

            //プーリー無しなのでP3は存在しない。
        }

        // P4,P5の座標計算
        if (this.orbitType == 1) {
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

            this.px[4] = -Math.sqrt(this.r0 ** 2 - this.d0 ** 2) + this.r3 * Math.cos(angle2);
            this.py[4] = -this.d0 - this.r3 * Math.sin(angle2);

            this.px[5] = this.r2 * Math.cos(angle2) - this.dlcs;
            this.py[5] = this.r2 * Math.sin(angle2);

        } else if (this.orbitType == 2) {

            //プーリー無しなのでP4は存在しない。

            this.px[5] = this.px[6];
            this.py[5] = -this.py[6];
        }
        
        //P1→P2
        this.x[1] = this.px[1];
        this.y[1] = this.py[1];
        setOnCircumference(0, 0, this.r1, this.px[1], this.py[1], this.px[2], this.py[2], this.x, this.y, this.cp, true);

        // P2→P5
        if (this.orbitType == 1) {
            //P2→P3
            tangent(this.px[2], this.py[2], this.px[3], this.py[3], this.x, this.y, this.cp);

            //P3→P4
            setOnCircumference(
            -Math.sqrt(this.r0 ** 2 - this.d0 ** 2), -this.d0, this.r3, this.px[3], this.py[3], this.px[4], this.py[4], this.x, this.y, this.cp, false);

            //P4→P5
            tangent(this.px[4], this.py[4], this.px[5], this.py[5], this.x, this.y, this.cp);

        } else if (this.orbitType == 2) {
            // P2→P5
            tangent(this.px[2], this.py[2], this.px[5], this.py[5], this.x, this.y, this.cp);
        }
        
        //P5→P6
        setOnCircumference(-this.dlcs, 0, this.r2, this.px[5], this.py[5], this.px[6], this.py[6], this.x, this.y, this.cp, true);

        //P6→P1
        tangent(this.px[6], this.py[6], this.px[1], this.py[1], this.x, this.y, this.cp);

        //x[0],y[0]に最終コマの座標をセット
        this.x[0] = this.x[this.x.length - 1];
        this.y[0] = this.y[this.y.length - 1];


        if (this.orbitType == 2) {
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
        }
    }

    //最長チェーンステー長計算_Sim2
    setOrbit_dlcsMax(nc2) {
    
        // チェーンステー長導出
        let dlt = 1.0;
    
        for (let nnn = 0; nnn <= 10000; nnn++) {
            this.dlcs = this.dlcs + dlt;
        
            this.px = [];
            this.py = [];
            this.x = [];
            this.y = [];
        
            for (let mmm = 0; mmm <= 6; mmm++) {
                this.px[mmm] = 0;
                this.py[mmm] = 0;
            }
        
            this.orbitType = 2;
            this.setOrbit();
            let nc = this.x.length - 1;
        
            if (nc == nc2 + 1) {
                if (dlt == 1.0) {
                    this.dlcs = this.dlcs - 4*dlt;
                    dlt = 0.1;
                } else if (dlt == 0.1) {
                    this.dlcs = this.dlcs - 4*dlt;
                    dlt = 0.05;
                } else if (dlt == 0.05) {
                    this.dlcs = this.dlcs - 4*dlt;
                    dlt = 0.01;
                } else if ((dlt == 0.01) && (nc2+1 == nc)) {
                    this.dlcs = this.dlcs - 4*dlt;
                    dlt = 0.001;
                }

            } else if ((dlt == 0.001) && (nc2 == nc)
                && (Math.abs((this.x[1] - this.x[0]) ** 2 + (this.y[1] - this.y[0]) ** 2 - this.cp ** 2) < 0.1)) {
                    break;
            } else if (nc > nc2 + 1) {
                return false;
            }
        }
        return true;
    }

    // 中心間最小距離(チェーンステー⇔チェーンピン)最適化_Sim1
    set_d0Min(nc2) {

        let nc;
    
        this.orbitType = 1;
        this.setOrbit();
    
        nc = this.x.length - 1;
    
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
            
            this.orbitType = 1;
            this.setOrbit();
            nc = this.x.length - 1;
        
            if (nc <= nc2) {
                if (dlt == 10.0) {
                    this.d0 = this.d0 - dlt;
                    dlt = 1.0;
                } else if (dlt == 1.0) {
                    this.d0 = this.d0 - dlt;
                    dlt = 0.1;
                } else if (dlt == 0.1) {
                    this.d0 = this.d0 - dlt;
                    dlt = 0.01;
                } else if (dlt == 0.01) {
                    this.d0 = this.d0 - dlt;
                    dlt = 0.001;
                } else if (dlt == 0.001) {
                    if (nc < nc2) {
                        return false;	//最適化失敗
                    } else if (Math.abs((this.x[1] - this.x[0]) ** 2 + (this.y[1] - this.y[0]) ** 2 - this.cp ** 2) < 0.1) {
                        break;
                    }
                } 
            }
        }
        return true;
    }

    // チェーンコマ数の最適化
    optimumLinkNumber(){
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
            nc = nc - 1; //標準ピッチでコマ数が奇数の時は1コマ詰める
        }
        }
        return nc;
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

//共通関数

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

    if (varFormat(dist1 + dist2, "0.000") == varFormat(dist0, "0.000")) {
        return true;
    } else {
        return false;
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