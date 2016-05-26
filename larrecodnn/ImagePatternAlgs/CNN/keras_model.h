#ifndef KERAS_MODEL__H
#define KERAS_MODEL__H

#include <string>
#include <vector>
#include <fstream>
#include <iostream>

class DataChunk {
public:
  virtual ~DataChunk() {}
  virtual std::vector<float> get_1d() {
    return std::vector<float>();
  };
  virtual std::vector<std::vector<std::vector<float> > > get_3d() {
    return   std::vector<std::vector<std::vector<float> > > ();
  };

  virtual void set_data(std::vector<std::vector<std::vector<float> > > ) {};
  virtual void set_data(std::vector<float> ) {};
  //virtual unsigned int get_count();
  virtual void read_from_file(const std::string &fname) {};
  virtual void show_name() = 0;
  virtual void show_values() = 0;
};

class DataChunk2D : public DataChunk {
public:
  std::vector<std::vector<std::vector<float> > > get_3d() {
    return data;
  };
  virtual void set_data(std::vector<std::vector<std::vector<float> > > d) { data = d; };

  void show_name() {
    std::cout << "DataChunk2D " << data.size() << "x" << data[0].size() << "x" << data[0][0].size() << std::endl;
  }

  void show_values() {
    std::cout << "DataChunk2D values:" << std::endl;
    for(size_t i = 0; i < data.size(); ++i) {
      std::cout << "Kernel " << i << std::endl;
      for(size_t j = 0; j < data[0].size(); ++j) {
        for(size_t k = 0; k < data[0][0].size(); ++k) {
          std::cout << data[i][j][k] << " ";
        }
        std::cout << std::endl;
      }
    }
  }
  //unsigned int get_count() {
  //  return data.size()*data[0].size()*data[0][0].size();
  //}

  void read_from_file(const std::string &fname);
  std::vector<std::vector<std::vector<float> > > data; // depth, rows, cols

  int m_depth;
  int m_rows;
  int m_cols;
};

class DataChunkFlat : public DataChunk {
public:
  void set_data(std::vector<float> d) { f = d; };

  void show_name() {
    std::cout << "DataChunkFlat " << f.size() << std::endl;
  }
  void show_values() {
    std::cout << "DataChunkFlat values:" << std::endl;
    for(size_t i = 0; i < f.size(); ++i) std::cout << f[i] << " ";
    std::cout << std::endl;
  }
  void read_from_file(const std::string &fname) {};
  std::vector<float> f;
  std::vector<float> get_1d() {
    return f;
  };
  //unsigned int get_count() { return f.size(); }
};

class Layer {
public:
  virtual void load_weights(std::ifstream &fin) = 0;
  virtual DataChunk* compute_output(DataChunk*) = 0;

  Layer(std::string name) : m_name(name) {}
  virtual ~Layer() {}

  virtual unsigned int get_input_rows() const = 0;
  virtual unsigned int get_input_cols() const = 0;

  std::string get_name() { return m_name; }
  std::string m_name;
};


class LayerFlatten : public Layer {
public:
  LayerFlatten() : Layer("Flatten") {}
  void load_weights(std::ifstream &fin) {};
  DataChunk* compute_output(DataChunk*);

  virtual unsigned int get_input_rows() const { return 0; } // !!!!! remember to implement !!!!
  virtual unsigned int get_input_cols() const { return 0; } // !!!!! remember to implement !!!!
};


class LayerMaxPooling : public Layer {
public:
  LayerMaxPooling() : Layer("MaxPooling2D") {};

  void load_weights(std::ifstream &fin);
  DataChunk* compute_output(DataChunk*);

  virtual unsigned int get_input_rows() const { return 0; } // !!!!! remember to implement !!!!
  virtual unsigned int get_input_cols() const { return 0; } // !!!!! remember to implement !!!!

  int m_pool_x;
  int m_pool_y;

};

class LayerActivation : public Layer {
public:
  LayerActivation() : Layer("Activation") {}
  void load_weights(std::ifstream &fin);
  DataChunk* compute_output(DataChunk*);

  virtual unsigned int get_input_rows() const { return 0; } // !!!!! remember to implement !!!!
  virtual unsigned int get_input_cols() const { return 0; } // !!!!! remember to implement !!!!

  std::string m_activation_type;
};

class LayerConv2D : public Layer {
public:
  LayerConv2D() : Layer("Conv2D") {}

  void load_weights(std::ifstream &fin);
  DataChunk* compute_output(DataChunk*);
  std::vector<std::vector<std::vector<std::vector<float> > > > m_kernels; // kernel, depth, rows, cols
  std::vector<float> m_bias; // kernel

  virtual unsigned int get_input_rows() const { return m_rows; }
  virtual unsigned int get_input_cols() const { return m_cols; }

  int m_kernels_cnt;
  int m_depth;
  int m_rows;
  int m_cols;
};

class LayerDense : public Layer {
public:
  LayerDense() : Layer("Dense") {}

  void load_weights(std::ifstream &fin);
  DataChunk* compute_output(DataChunk*);
  std::vector<std::vector<float> > m_weights; //input, neuron
  std::vector<float> m_bias; // neuron

  virtual unsigned int get_input_rows() const { return 1; } // flat, just one row
  virtual unsigned int get_input_cols() const { return m_input_cnt; }

  int m_input_cnt;
  int m_neurons;
};

class KerasModel {
public:
  KerasModel(const std::string &input_fname);
  ~KerasModel();
  std::vector<float> compute_output(DataChunk *dc);

  unsigned int get_input_rows() const { return m_layers.front()->get_input_rows(); }
  unsigned int get_input_cols() const { return m_layers.front()->get_input_cols(); }
  int get_output_length() const { return 0; } // !!!!! remember to implement !!!!

private:

  void load_weights(const std::string &input_fname);
  int m_layers_cnt; // number of layers
  std::vector<Layer *> m_layers; // container with layers

};

#endif
