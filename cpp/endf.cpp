#ifndef PYNE_IS_AMALGAMATED
  #include "endf_base.h"
  #include "endf.h"
#endif

bool pyne::endf::operator <(const pyne::endf::endf_id &lhs, const pyne::endf::endf_id &rhs) {
  if (lhs.mf < rhs.mf)
    return true;
  else if(lhs.mf > rhs.mf)
    return false;
  else if(lhs.mt < rhs.mt)
    return true;
  else if(lhs.mt > rhs.mt)
    return false;
  else if(lhs.mat < rhs.mat)
    return true;
  else
    return false;
}

bool pyne::endf::operator ==(const pyne::endf::endf_id &lhs, const pyne::endf::endf_id &rhs) {
  if ((lhs.mf == rhs.mf) && (lhs.mt == rhs.mt) && (lhs.mat == rhs.mat))
    return true;
  else
    return false;
}

pyne::endf::library::~library() {
  std::map<endf_id, mt_base*>::iterator endf_begin, endf_end, it;
  endf_begin = contents.begin();
  endf_end = contents.end();
  for (it = endf_begin; it!= endf_end; ++it) {
      delete it->second;
  }
}

pyne::endf::endf_id pyne::endf::make_endf_id(mt_base input) {
  endf_id ret;
  ret.mat = input.mat;
  ret.mf = input.mf;
  ret.mt = input.mt;
  return ret;
}

std::vector<std::vector<int> >  pyne::endf::library::get_content_list() {
  std::vector<std::vector<int> > retlist;
  std::map<endf_id, mt_base*>::iterator endf_begin, endf_end, it;
  endf_begin = contents.begin();
  endf_end = contents.end();
  for (it = endf_begin; it!= endf_end; ++it) {
      std::vector<int> id;
      id.push_back(it->first.mat);
      id.push_back(it->first.mf);
      id.push_back(it->first.mt);
      retlist.push_back(id);
  }
  return retlist;
}

void pyne::endf::library::read_endf(std::string filenm) {
  std::ifstream infile;
  std::string line;
  infile.open(filenm.c_str());
  getline(infile, line);// revision and other stuff

  int mf, mt, count, place;
  place = infile.tellg();
  getline(infile, line);
  infile.seekg(place);
  mt = atoi(line.substr(72, 3).c_str());
  mt_451 *mt451;
  if (mt == 451){
    mt451 = new mt_451;
    *mt451 = read_451(infile);
    endf_id id = make_endf_id((mt_base) *mt451);
    contents[(endf_id) id] = mt451;
  }else
    return;

  std::vector<int> plist;
  plist.push_back(81);
  infile.seekg(place);
  for (int i = 0; i < mt451->nxc; ++i) {
    place = infile.tellg();
    getline(infile, line);
    infile.seekg(place);
    mt = atoi(line.substr(72, 3).c_str());
    count = 0;
    while (mt == 0) {
      count++;
      place = infile.tellg();
      getline(infile, line);
      mt = atoi(line.substr(72, 3).c_str());
    }
    if (count > 0)
      --count;
    plist.back() = plist.back() + count*81;
    infile.seekg(place);
    int last = plist.back();
    plist.push_back(last + mt451->mt_list[i][2]*81 + 81);
    infile.seekg(plist.back());
  }
  for (int i = 1; i < plist.size() - 1; ++i){
    infile.clear();// reset any flags if a parser went nuts
    infile.seekg(plist[i]);
    mf = mt451->mt_list[i][0];
    mt = mt451->mt_list[i][1];
    if (mf == 1) {
      if (mt == 452) {
        mt_452_1 *temp = new mt_452_1;
        *temp = read_452_1(infile);
        endf_id id = make_endf_id((mt_base) *temp);
        contents[(endf_id) id] = temp;
      }else if (mt == 455) {
        mt_455_1 *temp = new mt_455_1;
        *temp = read_455_1(infile);
        endf_id id = make_endf_id((mt_base) *temp);
        contents[(endf_id) id] = temp;
      }else if (mt == 456) {
        mt_456_1 *temp = new mt_456_1;
        *temp = read_456_1(infile);
        endf_id id = make_endf_id((mt_base) *temp);
        contents[(endf_id) id] = temp;
      }else if (mt == 458) {
        mt_458_1 *temp = new mt_458_1;
        *temp = read_458_1(infile);
        endf_id id = make_endf_id((mt_base) *temp);
        contents[(endf_id) id] = temp;
      }else if (mt == 460) {
        mt_460_1 *temp = new mt_460_1;
        *temp = read_460_1(infile);
        endf_id id = make_endf_id((mt_base) *temp);
        contents[(endf_id) id] = temp;
      }
    }else if (mf == 8) {
      if ((mt == 454) || (mt == 459)) {
        mt_fpy_8 *temp = new mt_fpy_8;
        *temp = read_fpy_8(infile);
        endf_id id = make_endf_id((mt_base) *temp);
        contents[(endf_id) id] = temp;
      }
    }
  }
}

template<typename T> T pyne::endf::library::get (endf_id comp){
  std::map<endf_id, mt_base*>::iterator endf_iter, endf_end;
  endf_iter = contents.find(comp);
  endf_end = contents.end();
  if (endf_iter != endf_end) {
    return *dynamic_cast<T *>(endf_iter->second);
  }
  T ret;
  return ret;
};

template pyne::endf::mt_451 pyne::endf::library::get(pyne::endf::endf_id comp);
template pyne::endf::mt_452_1 pyne::endf::library::get(pyne::endf::endf_id comp);
template pyne::endf::mt_455_1 pyne::endf::library::get(pyne::endf::endf_id comp);
template pyne::endf::mt_456_1 pyne::endf::library::get(pyne::endf::endf_id comp);
template pyne::endf::mt_458_1 pyne::endf::library::get(pyne::endf::endf_id comp);
template pyne::endf::mt_460_1 pyne::endf::library::get(pyne::endf::endf_id comp);
template pyne::endf::mt_fpy_8 pyne::endf::library::get(pyne::endf::endf_id comp);

template<typename T> T pyne::endf::library::get (int mat, int mf, int mt){
  endf_id tmp;
  tmp.mat = mat;
  tmp.mf = mf;
  tmp.mt = mt;
  return get<T>(tmp);
};

template pyne::endf::mt_451 pyne::endf::library::get(int mat, int mf, int mt);
template pyne::endf::mt_452_1 pyne::endf::library::get(int mat, int mf, int mt);
template pyne::endf::mt_455_1 pyne::endf::library::get(int mat, int mf, int mt);
template pyne::endf::mt_456_1 pyne::endf::library::get(int mat, int mf, int mt);
template pyne::endf::mt_458_1 pyne::endf::library::get(int mat, int mf, int mt);
template pyne::endf::mt_460_1 pyne::endf::library::get(int mat, int mf, int mt);
template pyne::endf::mt_fpy_8 pyne::endf::library::get(int mat, int mf, int mt);

template<typename T> std::vector<T> pyne::endf::library::getl (int mf, int mt){
  std::map<endf_id, mt_base*>::iterator endf_lower, endf_upper, it;
  endf_id tmp1, tmp2;
  tmp1.mat = 0;
  tmp2.mat = 9999999;
  tmp1.mf = mf;
  tmp1.mt = mt;
  tmp2.mf = mf;
  tmp2.mt = mt;
  endf_lower = contents.lower_bound(tmp1);
  endf_upper = contents.upper_bound(tmp2);
  std::vector<T> result;
  T ret;
  for (it = endf_lower; it!= endf_upper; ++it) {
    ret = *dynamic_cast<T *>(it->second);
    result.push_back(ret);
  }
  return result;
};


template std::vector<pyne::endf::mt_451> pyne::endf::library::getl(int mf, int mt);
template std::vector<pyne::endf::mt_452_1> pyne::endf::library::getl(int mf, int mt);
template std::vector<pyne::endf::mt_455_1> pyne::endf::library::getl(int mf, int mt);
template std::vector<pyne::endf::mt_456_1> pyne::endf::library::getl(int mf, int mt);
template std::vector<pyne::endf::mt_458_1> pyne::endf::library::getl(int mf, int mt);
template std::vector<pyne::endf::mt_460_1> pyne::endf::library::getl(int mf, int mt);
template std::vector<pyne::endf::mt_fpy_8> pyne::endf::library::getl(int mf, int mt);
