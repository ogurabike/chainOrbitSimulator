const iterationCount=10000000;

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
    isError;

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
        let dlt;
        let delta;
        let deltamin;
        let deltamin0;
        let angle1;
        let angle2;

        let c1x;
        let c1y;
        let c2x;
        let c2y;
        let c3x;
        let c3y;
        let p2x;
        let p2y;
        let p3x;
        let p3y;
        let p4x;
        let p4y;
        let p5x;
        let p5y;
        
        deltamin0 = 100000000;

        this.isError = false;
        this.x = [];
        this.y = [];
        this.px = [];
        this.py = [];

        c1x = 0;
        c1y = 0;
        c2x = - this.dlcs;
        c2y = 0;
        c3x = - Math.sqrt(this.r0 ** 2 - this.d0 ** 2);
        c3y = - this.d0;

        // P6の座標計算
        this.px[6] = -(this.dlcs ** 2 + this.r2 * (this.r1 - this.r2)) / this.dlcs;
        this.py[6] = this.r2 * Math.sqrt(this.dlcs ** 2 - (this.r1 - this.r2) ** 2) / this.dlcs;
        
        // P1の座標計算
        this.px[1] = -this.r1 * (this.r1 - this.r2) / this.dlcs;
        this.py[1] = this.r1 * Math.sqrt(this.dlcs ** 2 - (this.r1 - this.r2) ** 2) / this.dlcs;

        // チェーンステー長を伸ばしてチェーンをパツパツにした時（プーリー無し軌道）のプーリーの位置を導出
        if (this.orbitType == 2) {
            dlt = 10.0;
            for (let nnn = 0; nnn <= iterationCount; nnn++) {
                
                if (this.pulleyPositionCheck(2) == false) {
                    this.d0 = this.d0 - 2.1 * dlt;
                    dlt = dlt/10;
                    if (Math.abs(varFormat(1/dlt,"0")) == iterationCount) {
                        this.d0 = this.d0 - dlt;
                        break;
                    }
                }

                this.d0 = this.d0 + dlt;

                if (nnn == iterationCount) {
                    console.log("Line:316 iterationCount max");
                    this.isError=true;
                    return;
                }
            }
        }

        // プーリー位置に矛盾がある場合は計算はしない
        if (this.orbitType == 1) {
            if (this.pulleyPositionCheck(1) == false) {
                this.isError=true;
                return;
            }
            if (this.pulleyPositionCheck(2) == false) {
                this.isError=true;
                return;
            }
        }

        // P2,P3の座標計算
        if (this.orbitType == 1) {
            dlt = 0.1;
            delta = 0;
            deltamin = deltamin0;
            //angle1 = Math.PI + Math.abs(Math.atan(this.d0 / Math.sqrt(this.r0 ** 2 - this.d0 ** 2)));
            angle1 = 1 / 2 * Math.PI

            //console.log(",mmm,dlt,angle1,delta,deltamin");

            for (let mmm = 1; mmm <= iterationCount; mmm++) {

                delta = Math.abs(
                    this.r3 
                    - Math.abs(
                        Math.sqrt(this.r0 ** 2 - this.d0 ** 2) * Math.cos(angle1) 
                        + this.d0 * Math.sin(angle1) 
                        + this.r1
                    )
                );

                //console.log(","+mmm+","+dlt+","+angle1+","+delta+","+deltamin);
                
                p2x = this.r1 * Math.cos(angle1);
                p2y = this.r1 * Math.sin(angle1);
                p3x = -Math.sqrt(this.r0 ** 2 - this.d0 ** 2) + this.r3 * Math.abs(Math.cos(angle1));
                p3y = -this.d0 + this.r3 * Math.abs(Math.sin(angle1));
                //C1を起点とした時のC3,P2,P3をつなぐベクトルの外積のZ座標の向きを評価する。
                //C1C3 X C1P2 → 正
                //C1C3 X C1P3 → 負
                if ((outerProduct(c1x,c1y,c3x,c3y,p2x,p2y) > 0)
                    && (outerProduct(c1x,c1y,c3x,c3y,p3x,p3y) < 0)){
                    if (delta <= deltamin) {
                        deltamin = delta;
                        angle1 = angle1 + dlt;
                    } else if (Math.abs(varFormat(1/dlt,"0")) == iterationCount) {
                        angle1 = angle1 - dlt;
                        //console.log(","+mmm+","+dlt+","+angle1+","+delta+","+deltamin);
                        break;
                    } else {
                        deltamin = delta;
                        dlt = - dlt/10;
                        angle1 = angle1 + dlt;
                    }
                } else {
                    angle1 = angle1 + dlt;
                }

                if (mmm == iterationCount) {
                    console.log("Line:141 iterationCount max / dlt = "+dlt+" / angle1 = "+angle1+" / delta = "+delta);
                    this.isError=true;
                    return;
                }
            }

            p2x = this.r1 * Math.cos(angle1);
            p2y = this.r1 * Math.sin(angle1);
            p3x = -Math.sqrt(this.r0 ** 2 - this.d0 ** 2) + this.r3 * Math.abs(Math.cos(angle1));
            p3y = -this.d0 + this.r3 * Math.abs(Math.sin(angle1));
            
            this.px[2] = p2x;
            this.py[2] = p2y;

            this.px[3] = p3x;
            this.py[3] = p3y;

            // console.log("P2,P3 deltamin="+deltamin);

        } else if (this.orbitType == 2) {
            this.px[2] = this.px[1];
            this.py[2] = -this.py[1];

            //プーリー無しなのでP3は存在しない。
        }

        // P4,P5の座標計算
        if (this.orbitType == 1) {

            if (this.d0 == this.r2 + this.r3) {
                angle2 = 3 / 2 * Math.PI;
            } else {
                dlt = 0.1;
                delta = 0;
                deltamin = deltamin0;

                if (this.d0 > this.r2 + this.r3){
                    angle2 = Math.atan(-this.py[6] / (this.px[6]+this.dlcs));
                } else {
                    angle2 = 3 / 2 * Math.PI + dlt;
                }

                //console.log(",mmm,dlt,angle2,delta,deltamin,pulleyPositionCheck");

                for (let mmm = 1; mmm <= iterationCount; mmm++) {
                    
                    p4x = - Math.sqrt(this.r0 ** 2 - this.d0 ** 2) + this.r3 * Math.cos(angle2);
                    p4y = - this.d0 - this.r3 * Math.sin(angle2);
                    p5x = this.r2 * Math.cos(angle2) - this.dlcs;
                    p5y = this.r2 * Math.sin(angle2);
                    
                    delta = Math.abs(
                        this.r3
                        - Math.abs(
                            (c3x + this.dlcs) * Math.cos(angle2)
                            + c3y * Math.sin(angle2)
                            - this.r2
                        )
                    );

                    //console.log(","+mmm+","+dlt+","+angle2+","+delta+","+deltamin+","+this.pulleyPositionCheck(2));

                    //C2起点とした時のC3,P4,P5をつなぐベクトルの外積のZ座標の向きを評価する。
                    //C2C3 X C2P4 → 正
                    //C2C3 X C2P5 → 負
                    if (outerProduct(c2x,c2y,c3x,c3y,p4x,p4y) > 0
                        && outerProduct(c2x,c2y,c3x,c3y,p5x,p5y) < 0) {
                            if (delta <= deltamin) {
                                deltamin = delta;
                                angle2 = angle2 + dlt;
                            } else if (Math.abs(varFormat(1/dlt,"0")) == iterationCount) {
                                angle2 = angle2 - dlt;
                                break;
                            } else {
                                deltamin = delta;
                                dlt = - dlt/10;
                                angle2 = angle2 + dlt;
                            }           
                    } else {
                        angle2 = angle2 + dlt;
                    }

                    if (mmm == iterationCount) {
                        console.log("Line:229 iterationCount max");
                        this.isError=true;
                        return;
                    }
                }
            }
            
            p4x = - Math.sqrt(this.r0 ** 2 - this.d0 ** 2) + this.r3 * Math.cos(angle2);
            p4y = - this.d0 - this.r3 * Math.sin(angle2);
            p5x = this.r2 * Math.cos(angle2) - this.dlcs;
            p5y = this.r2 * Math.sin(angle2);

            this.px[4] = p4x;
            this.py[4] = p4y;
            this.px[5] = p5x;
            this.py[5] = p5y;

        } else if (this.orbitType == 2) {

            //プーリー無しなのでP4は存在しない。

            this.px[5] = this.px[6];
            this.py[5] = -this.py[6];
        }
        
        //P1→P2
        this.x[1] = this.px[1];
        this.y[1] = this.py[1];
        setOnCircumference(
            0, 0, this.r1,
            this.px[1], this.py[1], this.px[2], this.py[2], this.x, this.y,
            this.cp, true);
        // if (this.orbitType == 1) {console.log("P1→P2終点:" + (this.x.length-1));}

        // P2→P5
        if (this.orbitType == 1) {
            //P2→P3
            tangent(this.px[2], this.py[2], this.px[3], this.py[3], this.x, this.y, this.cp);
            // console.log("P2→P3終点:" + (this.x.length-1));

            //P3→P4
            setOnCircumference(
                c3x, c3y, this.r3,
                this.px[3], this.py[3], this.px[4], this.py[4], this.x, this.y,
                this.cp, false);
            // console.log("P3→P4終点:" + (this.x.length-1));

            //P4→P5
            tangent(this.px[4], this.py[4], this.px[5], this.py[5], this.x, this.y, this.cp);
            // console.log("P4→P5終点:" + (this.x.length-1));

        } else if (this.orbitType == 2) {
            // P2→P5
            tangent(this.px[2], this.py[2], this.px[5], this.py[5], this.x, this.y, this.cp);
        }
        
        //P5→P6
        setOnCircumference(
            -this.dlcs, 0, this.r2,
            this.px[5], this.py[5], this.px[6], this.py[6], this.x, this.y,
            this.cp, true);
        // if (this.orbitType == 1) {console.log("P5→P6終点:" + (this.x.length-1));}

        //P6→P1
        tangent(this.px[6], this.py[6], this.px[1], this.py[1], this.x, this.y, this.cp);
        // if (this.orbitType == 1) {console.log("P6→P1終点:" + (this.x.length-1));}

        //x[0],y[0]に最終コマの座標をセット
        this.x[0] = this.x[this.x.length - 1];
        this.y[0] = this.y[this.y.length - 1];
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
    
        let dlt = - 0.1;
        let d= 0;
        let dmin=100000000000;

        let radmax = 3 / 2 * Math.PI
        let rad = radmax;
        this.d0 = - this.r0 * Math.sin(rad);

        for (let mmm = 1; mmm <= iterationCount; mmm++) {
            nc = 0;

            this.orbitType = 1;
            this.setOrbit();
            nc = this.x.length - 1;

            d = Math.sqrt(Math.abs(
                (this.x[1] - this.x[0]) ** 2 
                + (this.y[1] - this.y[0]) ** 2
                - this.cp ** 2));

            if (!this.pulleyPositionCheck(2)){
                rad = rad + dlt;
                this.d0 = - this.r0 * Math.sin(rad);
            } else {
                if ((dlt < 0 && nc < nc2) || (dlt > 0 && nc > nc2)) {
                    rad = rad + dlt;
                    this.d0 = - this.r0 * Math.sin(rad);
                } else if (nc == nc2) {
                    if (d <= dmin) {
                        dmin=d;
                        rad = rad + dlt;
                        this.d0 = - this.r0 * Math.sin(rad);
                    } else if (Math.abs(dlt) < (1/iterationCount)) {

                        rad = rad - dlt;
                        this.d0 = - this.r0 * Math.sin(rad);

                        this.orbitType = 1;
                        this.setOrbit();

                        return true; 
                    } else {
                        dmin = d;
                        dlt = - dlt/10.0;
                        rad = rad + dlt;
                        this.d0 = - this.r0 * Math.sin(rad);
                    }
                } else {
                    dmin = d;
                    dlt = - dlt/10.0;
                    rad = rad + dlt;
                    this.d0 = - this.r0 * Math.sin(rad);
                }
            } 
        }
        this.isError=true;
        return false;
    }
    
    //最長チェーンステー長計算_Sim2
    setOrbit_dlcsMax(nc2) {
    
        // チェーンステー長導出
        let dlt = 1.0;
        let d= 0;
        let dmin=100000000000;
    
        for (let nnn = 0; nnn <= iterationCount; nnn++) {
        
            this.orbitType = 2;
            this.setOrbit();
            let nc = this.x.length - 1;

            d = Math.abs(
                (this.x[1] - this.x[0]) ** 2
                + (this.y[1] - this.y[0]) ** 2- this.cp ** 2);
            
            if ((dlt > 0 && nc < nc2) || (dlt < 0 && nc > nc2)) {
                this.dlcs = this.dlcs + dlt;
            } else if (nc == nc2) {
                if (d <= dmin) {
                    dmin = d;
                    this.dlcs = this.dlcs + dlt;
                } else if (Math.abs(varFormat(1/dlt,"0")) == iterationCount) {
                    
                        this.orbitType = 2;
                        this.dlcs = this.dlcs - dlt;
                        this.setOrbit();
                        return true;
                } else {
                    dmin = d;
                    dlt = - dlt/10;
                    this.dlcs = this.dlcs + dlt;
                }
            } else {
                dmin = d;
                dlt = - dlt/10;
                this.dlcs = this.dlcs + dlt;
            }
        }
        this.isError=true;
        return false;
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
    //
    //■ chineLineNo=1 : 上側のチェーンの場合
    //d0を小さくしすぎて、接線P6P1をより上側にプーリー行ってしまうと、
    //上側チェーンとした側チェーンが干渉してしまう。
    //
    //■ chineLineNo=2 : 下側のチェーンの場合
    //d0を増やしすぎて、接線P6P1をチェーンステーの軸中心(y=0)で反転させた直線より
    //下側にプーリーが行ってしまうとチェーンとプーリーが接しなくなってしまう。
    pulleyPositionCheck(
        chineLineNo
    ) {
  
        let slope;
        let intercept;
        let p1x;
        let p1y;
        let p6x;
        let p6y;
        let c3x;
        let c3y;
    
        if ((this.r0 < this.d0) || (!(chineLineNo==1) && !(chineLineNo==2))) {
            return false;
        }

        p1x = -this.r1 * (this.r1 - this.r2) / this.dlcs;
        p1y = this.r1 * Math.sqrt(this.dlcs ** 2 - (this.r1 - this.r2) ** 2) / this.dlcs;
    
        p6x = -(this.dlcs ** 2 + this.r2 * (this.r1 - this.r2)) / this.dlcs;
        p6y = this.r2 * Math.sqrt(this.dlcs ** 2 - (this.r1 - this.r2) ** 2) / this.dlcs;
        
        if (chineLineNo==1) {
            slope = linearfunction(p1x, p1y, p6x, p6y,)["slope"];
            intercept = linearfunction(p1x, p1y, p6x, p6y,)["intercept"];
        } else if (chineLineNo==2) {
            slope = linearfunction(p1x, -p1y, p6x, -p6y,)["slope"];
            intercept = linearfunction(p1x, -p1y, p6x, -p6y,)["intercept"];
        }
        
    
        c3x = -Math.sqrt(this.r0 ** 2 - this.d0 ** 2);
        c3y = -this.d0;
    
        if ((slope * c3x - c3y + intercept) / Math.sqrt(slope ** 2 + 1) <= this.r3) {
            if (chineLineNo==1) {
                return false;   //上側は接していると異常
            } else if (chineLineNo==2) {
                return true;   //下側は接しているのが正常
            }
        } else {
            if (chineLineNo==1) {
                return true;   //上側は接していないのが正常
            } else if (chineLineNo==2) {
                return false;   //下側は接していないと異常
            }
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
    cp,
    isClockwise) {
        let s;
        let t;
        let u;
        let v;
        
        let rad;
        let angle;

        let delta;
        let deltamin;
        
        let rad1;
        let rad2;
        let rad12;
        let deltarad12;
        let ndeltarad12;

        let dlt;
        
        let nc;
        
        //円弧（円周）の両端の点から配置できる点の最大数を求める。
        rad1 = Math.acos(Math.abs(x1 - x0) / r);
        if (x1 >= x0) {
            if (y1 >= y0) {
                // 0≦rad1≦π/2
                rad1 = rad1;
            } else {
                // 3/2 π≦rad1＜2π
                rad1 = 2 * Math.PI - rad1;
            }
        } else {
            if (y1 >= y0) {
                // π/2＜rad1≦π
                rad1 = Math.PI - rad1;
            } else {
                // π＜rad1＜3/2 π
                rad1 = Math.PI + rad1;
            }
        }
        
        rad2 = Math.acos(Math.abs(x2 - x0) / r);
        if (x2 >= x0) {
            if (y2 >= y0) {
                // 0≦rad2≦π/2
                rad2 = rad2;
            } else {
                // 3/2 π≦rad2＜2π
                rad2 = 2 * Math.PI - rad2;
            }
        } else {
            if (y2 >= y0) {
                // π/2＜rad2≦π
                rad2 = Math.PI - rad2;
            } else {
                // π＜rad2＜3/2 π
                rad2 = Math.PI + rad2;
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

        // 底辺の長さcp、残りの二辺の長さrの二等辺三角形の頂角をΘとすると
        //　cp/2 = r sin(Θ/2) すなわち Θ = 2 asin (cp/2/r)
        //　チェーン一コマ毎に座標は中心(x0,y0)、半径rの円周上を上記Θ回転配置される。
        deltarad12 = Math.abs(2 * Math.asin(cp / 2 / r));
        ndeltarad12 = Math.floor(rad12 / deltarad12);
        
        //接線上の点(s,t)から指定距離cpの位置にある円周上の最初の点(u,v)を求める
        nc = x.length - 1;
        s = x[nc];
        t = y[nc];
        let deltamin0 = 10000000000;
        dlt = 0.1;
        delta = 0;
        deltamin = deltamin0;

        if (!(isClockwise)) {
            dlt = - dlt;
        }
        angle = rad1;

        for (let mmm = 1; mmm <= iterationCount; mmm++) {
            u = r * Math.cos(angle) + x0;
            v = r * Math.sin(angle) + y0;
            delta = Math.abs((s - u) ** 2 + (t - v) ** 2 - cp ** 2);

            if (mmm == iterationCount) {console.log("Line:518 iterationCount max");}
            
            if (delta <= deltamin) {
                deltamin = delta;
                angle = angle - dlt;
            } else if (Math.abs(varFormat(1/dlt,"0")) == iterationCount) {
                angle = angle + dlt;
                break;
            } else {
                angle = angle + 2 * dlt;
                dlt = dlt/10;
                u = r * Math.cos(angle) + x0;
                v = r * Math.sin(angle) + y0;
                deltamin = Math.abs((s - u) ** 2 + (t - v) ** 2 - cp ** 2);
            }
        }
        
        nc = nc + 1
        x[nc] = r * Math.cos(angle) + x0;
        y[nc] = r * Math.sin(angle) + y0;
        
        //残りの点を配置する。
        for (let mmm = 1; mmm <= ndeltarad12; mmm++) {
        
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

    if (!(isOnLineSegment(p1x, p1y, p2x, p2y, u, v))) {

        u = solv2(a, b, c)["x2"];
        v = slope * u + intercept;

        if (!(isOnLineSegment(p1x, p1y, p2x, p2y, u, v))) {
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

    //点(x1,y1)→点(x,y)の距離
    dist1 = Math.sqrt((x - x1) ** 2 + (y - y1) ** 2);

    //点(x2,y2)→点(x,y)の距離
    dist2 = Math.sqrt((x - x2) ** 2 + (y - y2) ** 2);

    if (varFormat(dist1 + dist2, "0.00000") == varFormat(dist0, "0.00000")) {
        return true;
    } else {
        return false;
    }

}

//外積計算 AB X AC Z座標を返す。
function outerProduct(ax,ay,bx,by,cx,cy) {
    return varFormat((bx-ax)*(cy-ay)-(cx-ax)*(by-ay),"0.00000")+0;   //最後+0しないと-0問題が出る
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