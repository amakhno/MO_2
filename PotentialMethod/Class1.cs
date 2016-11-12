﻿using System;
using System.Drawing;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PotentialMethod
{
    class FindWay
    {
        FindWay Father;
        Point Root;
        FindWay[] Childrens;
        Point[] mAllowed;
        Point Begining;
        //true - вниз/вверх
        //false - влево/вправо
        bool flag;
        public FindWay(int x, int y, bool _flag, Point[] _mAllowed, Point _Beg, FindWay _Father)
        {
            Begining = _Beg;
            flag = _flag;
            Root = new Point(x, y);
            mAllowed = _mAllowed;
            Father = _Father;
        }
        public Boolean BuildTree()
        {
            Point[] ps = new Point[mAllowed.Length];
            int Count = 0;
            for (int i = 0; i < mAllowed.Length; i++)
                if (flag)
                {
                    if (Root.Y == mAllowed[i].Y)
                    {
                        Count++;
                        ps[Count - 1] = mAllowed[i];
                    }

                }
                else
                    if (Root.X == mAllowed[i].X)
                {
                    Count++;
                    ps[Count - 1] = mAllowed[i];
                }

            FindWay fwu = this;
            Childrens = new FindWay[Count];
            //Point[] ss = new Point[mAllowed.Length];
            int k = 0;
            for (int i = 0; i < Count; i++)
            {
                if (ps[i] == Root) continue;
                if (ps[i] == Begining)
                {
                    while (fwu != null)
                    {
                        mAllowed[k] = fwu.Root;
                        fwu = fwu.Father;
                        k++;
                    };
                    for (; k < mAllowed.Length; k++) mAllowed[k] = new Point(-1, -1);
                    return true;
                }

                if (!Array.TrueForAll<Point>(ps, p => ((p.X == 0) && (p.Y == 0))))
                {
                    Childrens[i] = new FindWay(ps[i].X, ps[i].Y, !flag, mAllowed, Begining, this);
                    Boolean result = Childrens[i].BuildTree();
                    if (result) return true;
                }
            }
            return false;
        }

    }
}
