#pragma once
#include <QAbstractTableModel>
#include <vector>
#include <QString>

namespace cas::gui {

class MatrixTableModel : public QAbstractTableModel {
    Q_OBJECT
    Q_PROPERTY(int rows READ rowCount NOTIFY dimensionsChanged)
    Q_PROPERTY(int columns READ columnCount NOTIFY dimensionsChanged)

public:
    explicit MatrixTableModel(int rows = 3, int cols = 3, QObject* parent = nullptr)
        : QAbstractTableModel(parent), m_data(rows, std::vector<QString>(cols, "0")) {}

    int rowCount(const QModelIndex& parent = QModelIndex()) const override {
        return static_cast<int>(m_data.size());
    }

    int columnCount(const QModelIndex& parent = QModelIndex()) const override {
        return m_data.empty() ? 0 : static_cast<int>(m_data[0].size());
    }

    QVariant data(const QModelIndex& index, int role = Qt::DisplayRole) const override {
        if (!index.isValid() || role != Qt::DisplayRole) return QVariant();
        return m_data[index.row()][index.column()];
    }

    bool setData(const QModelIndex& index, const QVariant& value, int role = Qt::EditRole) override {
        if (index.isValid() && role == Qt::EditRole) {
            m_data[index.row()][index.column()] = value.toString();
            emit dataChanged(index, index, {Qt::DisplayRole, Qt::EditRole});
            return true;
        }
        return false;
    }

    Qt::ItemFlags flags(const QModelIndex& index) const override {
        if (!index.isValid()) return Qt::NoItemFlags;
        return Qt::ItemIsEditable | Qt::ItemIsEnabled | Qt::ItemIsSelectable;
    }

    Q_INVOKABLE void resize(int r, int c) {
        beginResetModel();
        m_data.assign(r, std::vector<QString>(c, "0"));
        endResetModel();
        emit dimensionsChanged();
    }

    Q_INVOKABLE QString toLatex() const {
        QString res = "\\begin{pmatrix} ";
        for (int i = 0; i < rowCount(); ++i) {
            for (int j = 0; j < columnCount(); ++j) {
                res += m_data[i][j];
                if (j < columnCount() - 1) res += " & ";
            }
            if (i < rowCount() - 1) res += " \\\\ ";
        }
        res += " \\end{pmatrix}";
        return res;
    }

signals:
    void dimensionsChanged();

private:
    std::vector<std::vector<QString>> m_data;
};

} // namespace cas::gui
