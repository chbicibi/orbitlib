#! ruby -Ku
# version 0.3

def load_satcat path
  # ファイルを開いて1行ずつ処理する
  # 戻り値はハッシュを要素に持つ配列
  # => [{:satellite_number => satellite_number, ...}, ...]
  # フォーマット参照: "https://celestrak.com/satcat/satcat-format.asp"
  File.foreach(path).map do |line|
    {
      satellite_number: line[13..17].to_i,
      satellite_name: line[23..46].strip,
      launch_date: line[56..65],
      orbital_period: line[87..93].to_f,
      inclination: line[96..100].to_f,
      apogee_altitude: line[103..108].to_f,
      perigee_altitude: line[111..116].to_f,
      rcs: line[119..126].to_f,
      data_type: :satcat,
      raw_txt: line
    }
  end
end

def load_tle path
  # ファイルを開いて3行ずつ処理する
  # 戻り値はハッシュを要素に持つ配列
  # => [{:satellite_name => satellite_name, ...}, ...]
  # フォーマット参照: "http://celestrak.com/columns/v04n03/"
  File.foreach(path).each_slice(3).map do |lines|
    {
      satellite_name: lines[0][0..23].strip,
      satellite_number: lines[1][2..6].to_i,
      epoch_year: lines[1][18..19].to_i,
      epoch_day: lines[1][20..31].to_f,
      inclination: lines[2][8..15].to_f,
      raan: lines[2][17..24].to_f,
      eccentricity: lines[2][26..32].insert(0, "0.").to_f,
      argument_of_perigee: lines[2][34..41].to_f,
      mean_anomaly: lines[2][43..50].to_f,
      mean_motion: lines[2][52..62].to_f,
      revolution_number_at_epoch: lines[2][63..67].to_i,
      data_type: :tle,
      raw_txt: lines.join
    }
  end
end

def classify_file path
  if File.open(path, "r"){ |file| file.readline }.size > 25
    :satcat
  elsif
    :tle
  else
    :unknown
  end
end

def merge_data list0, list1
  # 2つの配列を結合し，各要素の :satellite_number キーでグループ分けし，
  # 要素が2つのものを抽出し(tleとrcsの両方が揃っているものを残す)，各要素をハッシュに変換
  # 入力パラメータ順は任意
  # 戻り値は要素数2のハッシュを要素に持つ配列
  # => [{:tle => tle_data, :satcat => satcat_data}, ...]
  (list0 + list1).group_by{ |item| item[:satellite_number] }
                 .select{ |key, value| value.size == 2 }
                 .map do |key, values|
    values.each_with_object({}){ |item, hash| hash[item[:data_type]] = item }
  end
end

def output_satcat_data data_list, path:
  File.open(path, "w") do |file|
    data_list.each do |data|
      file.write data[:raw_txt]
    end
  end
end

def output_data data_list, tle_path:, rcs_path:
  # 2つのファイルを同時に開いて書き出す
  File.open(tle_path, "w") do |file_tle|
    File.open(rcs_path, "w") do |file_rcs|
      data_list.each do |data_pair|
        file_tle.write data_pair[:tle][:raw_txt]
        file_rcs.write "#{data_pair[:satcat][:rcs]}\n"
      end
    end
  end
end

def main
  # コマンドライン引数処理
  inputs, options = ARGV.partition{ |arg| File.file? arg }

  # 入出力ファイル名設定
  input_path_satcat = inputs.find{ |path| classify_file(path) == :satcat } || "satcat.txt"
  input_path_tle = inputs.find{ |path| classify_file(path) == :tle } || "tle_20171012.txt"
  output_path_rcs = "RCS_list.txt"
  output_path_tle = "debri_elements.txt"

  puts %Q(input(satcat) : "#{input_path_satcat}")
  puts %Q(input(tle)    : "#{input_path_tle}")
  puts %Q(output(satcat): "#{output_path_rcs}")
  puts %Q(output(tle)   : "#{output_path_tle}")

  raise "入力ファイルがありません" unless [input_path_satcat, input_path_tle].all?{ |path| File.file? path }

  pattern = options.find{ |arg| arg[/^name\K(.+)/] } && $1 || "FENGYUN 1C DEB"

  # SATCAT テキストファイルを読み込み名前が"FENGYUN 1C DEB"にマッチするものを抽出
  satcat_list = load_satcat(input_path_satcat).select{ |item| item[:satellite_name].match(/#{pattern}/) }

  # 2行軌道要素(TLE) テキストファイルを読み込み名前が"FENGYUN 1C DEB"にマッチするものを抽出
  tle_list = load_tle(input_path_tle).select{ |item| item[:satellite_name].match(/#{pattern}/) }

  # return output_satcat_data satcat_list, path: "satcat_FENGYUN_1C_DEB.txt"

  # SATCATデータ配列とTLEデータ配列をマージする
  data_list = merge_data satcat_list, tle_list

  # 出力データ選択
  data_list = ARGV.reduce(data_list) do |result, option|
    case option
    when "sort"
      # データ配列をRCS降順でソートする
      result.sort_by{ |item| -item[:satcat][:rcs] }
    when "reverse"
      # データ配列を逆順にする
      result.reverse
    when "shuffle"
      # データ配列をランダムに並び替える
      result.shuffle
    when /^take *(\d+)/
      # データ配列から先頭のn個を選択する
      result.take $1.to_i
    when /^sample *(\d+)/
      # データ配列からランダムにn個選択する(順序はバラバラになる)
      result.sample $1.to_i
    else
      result
    end
    # TODO: 不正な入力の際の処理を追加する
  end

  # ファイル書き出し
  output_data data_list, tle_path: output_path_tle, rcs_path: output_path_rcs
end

main

# 使い方(コマンドライン引数)
# 1. 引数なし => "satcat.txt", "tle.txt" ファイルを現在のディレクトリに入れる
# 2. 入力ファイル名をコマンドライン引数にする(入力が無かったものは1.のデフォルト値が使われる． 順序任意)
# 3. オプション引数:
#   3-1. sort => データ配列をRCS降順でソートする
#   3-2. reverse => データ配列を逆順にする
#   3-3. shuffle => データ配列をランダムに並び替える
#   3-4. take[n] => データ配列から先頭のn個を選択する(括弧，スペースは不要)
#   3-5. sample[n] => データ配列からランダムにn個選択する(括弧，スペースは不要)
#   3-6. name[pattern] => 衛星の名前が正規表現でpatternにマッチするものを選択する(括弧，スペースは不要)
#  オプション引数は入力順に適用される(nameは除く)

# 例1: extract => デフォルトのファイル名が使われる
# 例2: extract "ファイル1" "ファイル2" sample100 => ランダムに100個選択し，順序ランダムで出力
# 例3: extract "ファイル1" "ファイル2" sort take100 => RCS 降順に100個選択し，RCS昇順で出力
# 例4: extract "ファイル1" "ファイル2" sort take1000 sample100 sort reverse => RCS 降順に1000個選択し，その中から100個選択してRCS昇順で出力
# 例5: extract "ファイル1" "ファイル2" nameCOMMET => 衛星の名前がCOMMETを含むものを選択する

__END__
