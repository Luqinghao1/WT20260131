/*
 * 文件名: fittingchart.cpp
 * 文件作用: 拟合绘图管理类实现文件
 * 修改记录:
 * 1. [优化] initializeCharts: 双对数坐标系 (Log-Log) 初始化时，X轴和Y轴不再默认强制使用科学计数法 ("eb")，改为通用格式 ("gb")，并取消固定精度限制。
 * 2. [保持] 保留了之前的半对数图 (Horner) 修复和笛卡尔图设置。
 */

#include "fittingchart.h"
#include <cmath>
#include <algorithm>
#include <QDebug>

FittingChart::FittingChart(QObject *parent)
    : QObject(parent), m_plotLogLog(nullptr), m_plotSemiLog(nullptr), m_plotCartesian(nullptr), m_calculatedPi(0.0)
{
}

void FittingChart::initializeCharts(MouseZoom* logLog, MouseZoom* semiLog, MouseZoom* cartesian)
{
    m_plotLogLog = logLog;
    m_plotSemiLog = semiLog;
    m_plotCartesian = cartesian;

    // --- 1. 初始化双对数图 (Log-Log) ---
    if (m_plotLogLog) {
        m_plotLogLog->xAxis->setLabel("时间 Time (h)");
        m_plotLogLog->yAxis->setLabel("压差 & 导数 (MPa)");

        QSharedPointer<QCPAxisTickerLog> logTickerX(new QCPAxisTickerLog);
        logTickerX->setLogBase(10.0);
        m_plotLogLog->xAxis->setTicker(logTickerX);
        m_plotLogLog->xAxis->setScaleType(QCPAxis::stLogarithmic);

        // [修改] 取消默认科学计数法 ("eb")，改为通用格式 ("gb")
        m_plotLogLog->xAxis->setNumberFormat("gb");
        // [修改] 取消默认精度限制，让通用格式自动决定精度
        // m_plotLogLog->xAxis->setNumberPrecision(1);

        QSharedPointer<QCPAxisTickerLog> logTickerY(new QCPAxisTickerLog);
        logTickerY->setLogBase(10.0);
        m_plotLogLog->yAxis->setTicker(logTickerY);
        m_plotLogLog->yAxis->setScaleType(QCPAxis::stLogarithmic);

        // [修改] 取消默认科学计数法 ("eb")，改为通用格式 ("gb")
        m_plotLogLog->yAxis->setNumberFormat("gb");
        // [修改] 取消默认精度限制
        // m_plotLogLog->yAxis->setNumberPrecision(1);

        // 设置图例位置：左上
        if(m_plotLogLog->axisRect() && m_plotLogLog->axisRect()->insetLayout())
            m_plotLogLog->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop | Qt::AlignLeft);
    }

    // --- 2. 初始化笛卡尔图 (Cartesian) ---
    if (m_plotCartesian) {
        m_plotCartesian->xAxis->setLabel("时间 Time (h)");
        m_plotCartesian->yAxis->setLabel("压差 Delta P (MPa)");

        QSharedPointer<QCPAxisTicker> linearTicker(new QCPAxisTicker);
        m_plotCartesian->xAxis->setTicker(linearTicker);
        m_plotCartesian->xAxis->setScaleType(QCPAxis::stLinear);
        m_plotCartesian->xAxis->setNumberFormat("gb");

        m_plotCartesian->yAxis->setTicker(linearTicker);
        m_plotCartesian->yAxis->setScaleType(QCPAxis::stLinear);
        m_plotCartesian->yAxis->setNumberFormat("gb");

        // 设置图例位置：右下
        if(m_plotCartesian->axisRect() && m_plotCartesian->axisRect()->insetLayout())
            m_plotCartesian->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignBottom | Qt::AlignRight);
    }

    // --- 3. 半对数图 (Semi-Log / Horner) ---
    if (m_plotSemiLog) {
        // 设置图例位置：左上
        if(m_plotSemiLog->axisRect() && m_plotSemiLog->axisRect()->insetLayout())
            m_plotSemiLog->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop | Qt::AlignLeft);
    }
}

void FittingChart::setObservedData(const QVector<double>& t, const QVector<double>& deltaP,
                                   const QVector<double>& deriv, const QVector<double>& rawP)
{
    m_obsT = t;
    m_obsDeltaP = deltaP;
    m_obsDeriv = deriv;
    m_obsRawP = rawP;
}

void FittingChart::setSettings(const FittingDataSettings& settings)
{
    m_settings = settings;
}

void FittingChart::plotAll(const QVector<double>& t_model, const QVector<double>& p_model, const QVector<double>& d_model, bool isModelValid, bool autoScale)
{
    if (!m_plotLogLog || !m_plotSemiLog || !m_plotCartesian) return;

    plotLogLog(t_model, p_model, d_model, isModelValid, autoScale);
    plotSemiLog(t_model, p_model, d_model, isModelValid, autoScale);
    plotCartesian(t_model, p_model, d_model, isModelValid, autoScale);
}

void FittingChart::plotLogLog(const QVector<double>& tm, const QVector<double>& pm, const QVector<double>& dm, bool hasModel, bool autoScale)
{
    MouseZoom* plot = m_plotLogLog;
    plot->clearGraphs();
    plot->clearItems();

    // 实测数据
    QVector<double> vt, vp, vd;
    for(int i=0; i<m_obsT.size(); ++i) {
        if (m_obsT[i] > 1e-10 && m_obsDeltaP[i] > 1e-10) {
            vt << m_obsT[i];
            vp << m_obsDeltaP[i];
            if (i < m_obsDeriv.size() && m_obsDeriv[i] > 1e-10) vd << m_obsDeriv[i];
            else vd << 1e-10;
        }
    }

    plot->addGraph();
    plot->graph(0)->setData(vt, vp);
    plot->graph(0)->setPen(Qt::NoPen);
    plot->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QColor(0, 100, 0), 6));
    plot->graph(0)->setName("实测压差");

    plot->addGraph();
    plot->graph(1)->setData(vt, vd);
    plot->graph(1)->setPen(Qt::NoPen);
    plot->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssTriangle, Qt::magenta, 6));
    plot->graph(1)->setName("实测导数");

    // 理论曲线
    if (hasModel) {
        QVector<double> vtm, vpm, vdm;
        for(int i=0; i<tm.size(); ++i) {
            if(tm[i] > 1e-10) {
                vtm << tm[i];
                vpm << (pm[i] > 1e-10 ? pm[i] : 1e-10);
                vdm << (dm[i] > 1e-10 ? dm[i] : 1e-10);
            }
        }
        plot->addGraph();
        plot->graph(2)->setData(vtm, vpm);
        plot->graph(2)->setPen(QPen(Qt::red, 2));
        plot->graph(2)->setName("理论压差");

        plot->addGraph();
        plot->graph(3)->setData(vtm, vdm);
        plot->graph(3)->setPen(QPen(Qt::blue, 2));
        plot->graph(3)->setName("理论导数");
    }

    // 仅在 autoScale 为 true 时复位视图
    if (autoScale) {
        plot->rescaleAxes();
        plot->xAxis->scaleRange(1.1, plot->xAxis->range().center());
        plot->yAxis->scaleRange(1.1, plot->yAxis->range().center());
    }

    if (m_settings.testType == Test_Buildup && m_calculatedPi > 1e-6) {
        showResultOnLogPlot();
    }
}

void FittingChart::plotSemiLog(const QVector<double>& tm, const QVector<double>& pm, const QVector<double>& dm, bool hasModel, bool autoScale)
{
    MouseZoom* plot = m_plotSemiLog;
    plot->clearGraphs();
    Q_UNUSED(dm);

    double tp = m_settings.producingTime;

    // 判断是否使用 Horner 图逻辑 (通常生产时间 tp > 0 时使用)
    bool useHorner = (tp > 1e-5);

    // --- 1. 准备实测数据 (Observed Pressure) ---
    QVector<double> vt, vp;
    for(int i=0; i<m_obsT.size(); ++i) {
        double dt = m_obsT[i];
        if (dt < 1e-6) continue;

        // X轴计算
        double x_val;
        if (useHorner) {
            // Horner 时间比: (tp + dt) / dt
            double ratio = (tp + dt) / dt;
            if (ratio <= 0) continue;
            x_val = log10(ratio); // 取对数
        } else {
            // MDH / Drawdown 模式: 直接使用时间 (后续通过坐标轴 Log 属性控制)
            x_val = dt;
        }

        // Y轴数据：使用原始实测压力 (Raw Pressure)，修正之前的压差数据错误
        double y_val = 0.0;
        if (i < m_obsRawP.size()) y_val = m_obsRawP[i];

        vt << x_val;
        vp << y_val;
    }

    plot->addGraph();
    plot->graph(0)->setData(vt, vp);
    plot->graph(0)->setPen(Qt::NoPen);
    plot->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QColor(0, 0, 180), 5)); // 蓝色圆点表示压力
    plot->graph(0)->setName("实测压力");

    // --- 2. 准备理论曲线 (Model Pressure) ---
    // 模型计算返回的是压差 pm (Pressure Drop)
    // 需要根据试井类型将其转换为绝对压力 P(t)
    if (hasModel) {
        QVector<double> vtm, vpm;

        // 确定起始参考压力
        // 如果有实测数据，以第一个点的压力为基准，使曲线起点对齐，便于视觉对比
        // 否则使用设定的初始压力 Pi
        double p_start = 0.0;
        if (!m_obsRawP.isEmpty()) p_start = m_obsRawP.first();
        else p_start = m_settings.initialPressure;

        for(int i=0; i<tm.size(); ++i) {
            double dt = tm[i];
            if (dt < 1e-10) continue;

            double x_val;
            if (useHorner) {
                double ratio = (tp + dt) / dt;
                if (ratio <= 0) continue;
                x_val = log10(ratio);
            } else {
                x_val = dt;
            }

            double y_val;
            // 根据试井类型转换压差为压力
            if (m_settings.testType == Test_Drawdown) {
                // 降压试井 (Drawdown): 压力下降 P(t) = P_start - ΔP
                y_val = p_start - pm[i];
            } else {
                // 恢复试井 (Buildup): 压力上升 P(t) = P_start + ΔP
                y_val = p_start + pm[i];
            }

            vtm << x_val;
            vpm << y_val;
        }

        plot->addGraph();
        plot->graph(1)->setData(vtm, vpm);
        plot->graph(1)->setPen(QPen(Qt::red, 2));
        plot->graph(1)->setName("理论压力");
    }

    // --- 3. 绘制 Horner 拟合线 (仅在 Horner 模式下) ---
    // 计算并绘制直线段，用于推算 Pi
    if (useHorner) {
        m_calculatedPi = calculateHornerPressure();
        if (m_calculatedPi > 0 && !vt.isEmpty()) {
            double startX = vt.first(); // 数据左侧点（对应较大的时间比）
            // 绘制从数据点延伸到 log((tp+dt)/dt)=0 (即无限时间) 的直线
            QVector<double> lineX = {startX, 0.0};
            // 简单两点连线：左端点数据Y -> 右端点 Pi
            QVector<double> lineY = {(vt.size() > 0 ? vp.first() : m_calculatedPi), m_calculatedPi};

            plot->addGraph();
            plot->graph(plot->graphCount()-1)->setData(lineX, lineY);
            plot->graph(plot->graphCount()-1)->setPen(QPen(Qt::red, 2, Qt::DashLine));
            plot->graph(plot->graphCount()-1)->setName("Horner 拟合线");
        }
    }

    // --- 4. 配置坐标轴 ---
    if (useHorner) {
        plot->xAxis->setLabel("Horner 时间比 lg((tp+dt)/dt)");
        plot->yAxis->setLabel("实测压力 Pressure (MPa)");

        // Horner 图 X 轴使用线性刻度（因为我们已经手动取了 log10）
        QSharedPointer<QCPAxisTicker> linearTicker(new QCPAxisTicker);
        plot->xAxis->setTicker(linearTicker);
        plot->xAxis->setScaleType(QCPAxis::stLinear);
        // Horner 图习惯上 X 轴反向：左边是大数值（短时间），右边是 0（长时间）
        plot->xAxis->setRangeReversed(true);

        plot->yAxis->setTicker(linearTicker);
        plot->yAxis->setScaleType(QCPAxis::stLinear);
    } else {
        // MDH / 普通半对数图
        plot->xAxis->setLabel("时间 Time (h)");
        plot->yAxis->setLabel("实测压力 Pressure (MPa)");

        QSharedPointer<QCPAxisTickerLog> logTicker(new QCPAxisTickerLog);
        plot->xAxis->setTicker(logTicker);
        plot->xAxis->setScaleType(QCPAxis::stLogarithmic);
        plot->xAxis->setRangeReversed(false);

        QSharedPointer<QCPAxisTicker> linearTicker(new QCPAxisTicker);
        plot->yAxis->setTicker(linearTicker);
        plot->yAxis->setScaleType(QCPAxis::stLinear);
    }

    if (autoScale) {
        plot->rescaleAxes();
        // Horner 模式下，右边界通常锁定在 0
        if(useHorner) {
            plot->xAxis->setRangeUpper(plot->xAxis->range().upper);
            plot->xAxis->setRangeLower(0.0);
        }
    }
}

void FittingChart::plotCartesian(const QVector<double>& tm, const QVector<double>& pm, const QVector<double>& dm, bool hasModel, bool autoScale)
{
    MouseZoom* plot = m_plotCartesian;
    plot->clearGraphs();
    Q_UNUSED(dm);

    QVector<double> vt, vp;
    for(int i=0; i<m_obsT.size(); ++i) { vt << m_obsT[i]; vp << m_obsDeltaP[i]; }

    plot->addGraph();
    plot->graph(0)->setData(vt, vp);
    plot->graph(0)->setPen(Qt::NoPen);
    plot->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QColor(0, 100, 0), 6));
    plot->graph(0)->setName("实测压差");

    if (hasModel) {
        QVector<double> vtm, vpm;
        for(int i=0; i<tm.size(); ++i) { vtm << tm[i]; vpm << pm[i]; }
        plot->addGraph();
        plot->graph(1)->setData(vtm, vpm);
        plot->graph(1)->setPen(QPen(Qt::red, 2));
        plot->graph(1)->setName("理论压差");
    }

    if (autoScale) plot->rescaleAxes();
}

void FittingChart::plotSampledPoints(const QVector<double>& t, const QVector<double>& p, const QVector<double>& d)
{
    if (!m_plotLogLog) return;

    QCPGraph* gP = m_plotLogLog->addGraph();
    gP->setData(t, p);
    gP->setPen(Qt::NoPen);
    gP->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(QColor(0, 100, 0)), QBrush(QColor(0, 100, 0)), 6));
    gP->setName("抽样压差");

    QCPGraph* gD = m_plotLogLog->addGraph();
    gD->setData(t, d);
    gD->setPen(Qt::NoPen);
    gD->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssTriangle, QPen(Qt::magenta), QBrush(Qt::magenta), 6));
    gD->setName("抽样导数");
}

double FittingChart::calculateHornerPressure()
{
    if (m_obsT.isEmpty() || m_obsRawP.isEmpty() || m_settings.producingTime <= 0) return 0.0;
    QVector<double> X, Y;
    for(int i=0; i<m_obsT.size(); ++i) {
        double dt = m_obsT[i];
        if(dt > 1e-5 && i < m_obsRawP.size()) {
            double val = (m_settings.producingTime + dt) / dt;
            if(val > 0) { X << log10(val); Y << m_obsRawP[i]; }
        }
    }
    if (X.size() < 5) return 0.0;

    int nPoints = X.size();
    int fitCount = nPoints * 0.3;
    if (fitCount < 3) fitCount = nPoints;
    double sumX = 0, sumY = 0, sumXY = 0, sumXX = 0;
    int count = 0;
    for(int i = nPoints - fitCount; i < nPoints; ++i) {
        sumX += X[i]; sumY += Y[i]; sumXY += X[i] * Y[i]; sumXX += X[i] * X[i]; count++;
    }
    if (std::abs(count * sumXX - sumX * sumX) < 1e-9) return 0.0;
    double slope = (count * sumXY - sumX * sumY) / (count * sumXX - sumX * sumX);
    return (sumY - slope * sumX) / count;
}

double FittingChart::getCalculatedInitialPressure() const { return m_calculatedPi; }

void FittingChart::showResultOnLogPlot()
{
    if(!m_plotLogLog) return;
    QCPItemText *textLabel = new QCPItemText(m_plotLogLog);
    textLabel->setPositionAlignment(Qt::AlignTop|Qt::AlignRight);
    textLabel->position->setType(QCPItemPosition::ptAxisRectRatio);
    textLabel->position->setCoords(0.95, 0.05);
    textLabel->setText(QString("Horner推算Pi: %1 MPa").arg(m_calculatedPi, 0, 'f', 2));
    textLabel->setFont(QFont("Microsoft YaHei", 10, QFont::Bold));
    textLabel->setColor(Qt::red);
    textLabel->setBrush(QBrush(QColor(255, 255, 255, 200)));
    textLabel->setPadding(QMargins(5, 5, 5, 5));
    textLabel->setPen(QPen(Qt::black));
}
